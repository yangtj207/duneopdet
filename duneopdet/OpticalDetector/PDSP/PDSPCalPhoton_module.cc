////////////////////////////////////////////////////////////////////////
// Class:       PDSPCalPhoton
// Plugin Type: producer
// File:        PDSPCalPhoton_module.cc
//
// Generated at Tue Aug  9 15:32:19 2022 by Tingjun Yang using cetskelgen
// 
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/search_path.h"

// LArSoft includes
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "PDSMeanRMS.h"

#include "TH1D.h"
#include "TFile.h"
#include "TVirtualFFT.h"
#include "TF1.h"

#include <vector>
#include <memory>
#include <complex>
#include <fftw3.h>

double elecrespfun(double *x, double *par){
  double y = 0;
  if (x[0]<par[3]){
    y = 0;
  }
  else if (x[0]<par[4]){
    y = par[0]*(1-exp(-(x[0]-par[3])/par[1]));
  }
  else{
    y = par[0]*(1-exp(-(par[4]-par[3])/par[1]))*exp(-(x[0]-par[4])/par[2]);
  }
  return y;
}

namespace pdsp {
  class PDSPCalPhoton;
}

class pdsp::PDSPCalPhoton : public art::EDProducer {
public:
  explicit PDSPCalPhoton(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDSPCalPhoton(PDSPCalPhoton const&) = delete;
  PDSPCalPhoton(PDSPCalPhoton&&) = delete;
  PDSPCalPhoton& operator=(PDSPCalPhoton const&) = delete;
  PDSPCalPhoton& operator=(PDSPCalPhoton&&) = delete;

  // Required functions.
  void produce(art::Event& e) override;

private:

  art::InputTag fRawOpWaveformLabel;
  std::string   fResponseFile;
  std::string   fFiltFunc;
  std::vector<double> fPECalib;
  
  int ntbin;
  TH1D *wienerfilter;
  TH1D *elecres;
  std::vector<double> vwiener;
  std::vector<std::complex<double>> vinvres;
  TF1 *filtfun;

  std::unique_ptr<pds::PDSMeanRMS> pdsMeanRMS_;  //Tool to get pedestal and RMS
  float roi_threshold;

};

double CalcMedian(std::vector<short int> scores){

  size_t size = scores.size();

  if (size == 0){
    return 0;  // Undefined, really.
  }
  else{
    sort(scores.begin(), scores.end());
    if (size % 2 == 0){
      return (scores[size / 2 - 1] + scores[size / 2]) / 2;
    }
    else {
      return scores[size / 2];
    }
  }
}

pdsp::PDSPCalPhoton::PDSPCalPhoton(fhicl::ParameterSet const& p)
  : EDProducer{p},
    fRawOpWaveformLabel(p.get<art::InputTag>("RawOpWaveformLabel")),
    fResponseFile(p.get<std::string>("ResponseFile")),
    fFiltFunc(p.get<std::string>("FiltFunc")),
    fPECalib(p.get<std::vector<double>>("PECalib")),
    ntbin(2000),
    pdsMeanRMS_{art::make_tool<pds::PDSMeanRMS>(p.get<fhicl::ParameterSet>("pdsMeanRMS"))},
    roi_threshold(p.get<float>("roi_threshold"))
{
  // Wiener filtered waveform
  produces<std::vector<recob::OpWaveform>>("wiener");
  produces<art::Assns<raw::OpDetWaveform, recob::OpWaveform>>("wiener");
  // Apply Gaussian filter and remove electronics response
  produces<std::vector<recob::OpWaveform>>();
  produces<art::Assns<raw::OpDetWaveform, recob::OpWaveform>>();

  // Get response functions
  std::string fullname;
  cet::search_path sp("FW_SEARCH_PATH");
  sp.find_file(fResponseFile, fullname);
  TFile *resfile;
  if (fullname.empty()) {
    throw cet::exception("PDSPCalPhoton") << "Unable to find response file "  << fResponseFile;
  }
  else{
    resfile = TFile::Open(fullname.c_str());
    wienerfilter = (TH1D*)resfile->Get("wienerfilter");
    elecres = (TH1D*)resfile->Get("elecres");
    if (!wienerfilter) throw cet::exception("PDSPCalPhoton") << "Unable to find wienerfilter";
    if (!elecres) throw cet::exception("PDSPCalPhoton") << "Unable to find elecres";
  }

  // Save wiener filter
  for (int i = 0; i<wienerfilter->GetNbinsX(); ++i){
    vwiener.push_back(wienerfilter->GetBinContent(i+1));
  }

  // Save inverse response function
  TF1 *respfun = new TF1("respfun", elecrespfun, 0, ntbin, 5);
  respfun->SetParameter(0, 10.23);
  respfun->SetParameter(1, 8.101);
  respfun->SetParameter(2, 66.43);
  respfun->SetParameter(3, 72.98);
  respfun->SetParameter(4, 94.63);

  fftw_complex *in, *out;
  fftw_plan pf;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  pf = fftw_plan_dft_1d(ntbin, in, out, FFTW_FORWARD, FFTW_MEASURE);
  double resp_integral = 0;
  for (int i = 0; i<ntbin; ++i){
    in[i][0] = respfun->Eval(i);
    in[i][1] = 0;
    resp_integral += in[i][0];
  }
  for (int i = 0; i<ntbin; ++i){
    in[i][0] *= resp_integral;
  }
  fftw_execute(pf);
  for (int i = 0; i<ntbin; ++i){
    std::complex<double> tmp (1,0);
    std::complex<double> res (out[i][0], out[i][1]);
    vinvres.push_back(tmp/res);
    //std::cout<<i<<" "<<real(vinvres.back())<<" "<<imag(vinvres.back())<<std::endl;
  }
  fftw_destroy_plan(pf);
  fftw_free(in); fftw_free(out);
  delete respfun;

  // Construct filter function
  filtfun = new TF1("filtfun",fFiltFunc.c_str(),0,150);

//  for (int i = 0; i<300; ++i){
//    std::cout<<i<<" "<<fPECalib[i]<<std::endl;
//  }
}

void pdsp::PDSPCalPhoton::produce(art::Event& e)
{
  // Implementation of required member function here.
  std::vector < art::Ptr < raw::OpDetWaveform > > wfList;
  auto wfHandle = e.getHandle<std::vector<raw::OpDetWaveform>>(fRawOpWaveformLabel);
  if (!wfHandle.isValid()){
    throw cet::exception("PDSPCalPhoton") << "Unable to get waveforms using label "  << fRawOpWaveformLabel;
  }
  else{
    art::fill_ptr_vector(wfList, wfHandle);
  }

  auto out_recowaveforms1 = std::make_unique< std::vector< recob::OpWaveform > >();
  auto out_recowaveforms2 = std::make_unique< std::vector< recob::OpWaveform > >();
  auto recorawassn1 = std::make_unique< art::Assns<raw::OpDetWaveform, recob::OpWaveform> >();
  auto recorawassn2 = std::make_unique< art::Assns<raw::OpDetWaveform, recob::OpWaveform> >();

  fftw_complex *in, *out;
  fftw_plan pf;
  fftw_plan pb;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntbin);
  pf = fftw_plan_dft_1d(ntbin, in, out, FFTW_FORWARD, FFTW_MEASURE);
  pb = fftw_plan_dft_1d(ntbin, in, out, FFTW_BACKWARD, FFTW_MEASURE);

  // Vector after applying Wiener filter
  std::vector<float> wfwiener(ntbin);
  // Vector after applying Gaussian filter and remove electronics response
  std::vector<float> wfimpulse(ntbin);

  for (auto const& wf : wfList) {
    if (int(wf->Waveform().size()) != ntbin){
      throw cet::exception("PDSPCalPhoton") << "Waveform size does not match "  << wf->Waveform().size()<<" "<<ntbin;
    }
    
    // Get baseline
    pdsMeanRMS_->Analyze(wf->Waveform());
    double baseline = pdsMeanRMS_->GetMean();

    for (int i = 0; i<ntbin; ++i){
      in[i][0] = wf->Waveform()[i] - baseline;
      in[i][1] = 0;
    }

    fftw_execute(pf);

    // FFT of waveform
    std::vector<std::complex<double>> wffft(wf->Waveform().size());

    for (int i = 0; i<ntbin; ++i){
      double filt = vwiener[i*vwiener.size()/ntbin];
      wffft[i] = std::complex<double>(out[i][0], out[i][1]);
      // Apply Wiener filter
      in[i][0] = out[i][0]*filt;
      in[i][1] = out[i][1]*filt;
    }

    fftw_execute(pb);

    for (int i = 0; i<ntbin; ++i){
      wfwiener[i] = out[i][0]/ntbin;
    }
    recob::OpWaveform out_recowaveFinal(wf->TimeStamp(), wf->ChannelNumber(), wfwiener);
    out_recowaveforms1->emplace_back(std::move(out_recowaveFinal));
    util::CreateAssn(*this, e, *out_recowaveforms1, wf, *recorawassn1, "wiener");

    // Apply Gaussian filter and remove electronics response
    for (int i = 0; i<ntbin; ++i){
      double filt = filtfun->Eval(150.*i/ntbin);
      if (i>ntbin/2){
        filt = filtfun->Eval(150.*(ntbin-i)/ntbin);
      }
      std::complex<double> tmp = wffft[i]*vinvres[i]*fPECalib[wf->ChannelNumber()];
      in[i][0] = real(tmp)*filt;
      in[i][1] = imag(tmp)*filt;
      //std::cout<<i<<" "<<wffft[i]<<" "<<vinvres[i*vinvres.size()/ntbin]<<" "<<tmp<<" "<<filt<<" "<<in[i][0]<<" "<<in[i][1]<<std::endl;
    }
    fftw_execute(pb);

    for (int i = 0; i<ntbin; ++i){
      wfimpulse[i] = out[i][0]/ntbin;
      //if (wf->ChannelNumber()==101) std::cout<<out[i][0]<<std::endl;
    }
    std::vector<float> wfimpulse_shift(ntbin);
    for (int i = 0; i<ntbin; ++i){
      if (i<95){
        wfimpulse_shift[i] = wfimpulse[ntbin-95+i];
      }
      else{
        wfimpulse_shift[i] = wfimpulse[i-95]; 
      }
    }
    std::vector<bool> inroi(ntbin,false);
    for (int i = 0; i<ntbin; ++i){
      if (wfimpulse_shift[i] > roi_threshold){
        for (int j = i - 100; j < i + 100; ++j){
          if (j>=0 && j<ntbin){
            inroi[j] = true;
          }
        }
      }
    }

    recob::OpWaveform::RegionsOfInterest_t rois(ntbin);
    std::vector<float> sigs;
    int lastsignaltick = -1;
    int roistart = -1;
    bool hasROI = false;
    for (int i = 0; i < ntbin; ++i) {
      if (inroi[i]) {
        hasROI = true;
        if (sigs.empty()) {
          sigs.push_back(wfimpulse_shift[i]);
          lastsignaltick = i;
          roistart = i;
        }
        else {
          if (i != lastsignaltick + 1) {
            rois.add_range(roistart, std::move(sigs));
            sigs.clear();
            sigs.push_back(wfimpulse_shift[i]);
            lastsignaltick = i;
            roistart = i;
          }
          else {
            sigs.push_back(wfimpulse_shift[i]);
            lastsignaltick = i;
          }
        }
      }
    }    
    if (!sigs.empty()) { rois.add_range(roistart, std::move(sigs)); }


    //recob::OpWaveform out_recowaveFinal2(wf->TimeStamp(), wf->ChannelNumber(), wfimpulse_shift);
    //out_recowaveforms2->emplace_back(std::move(out_recowaveFinal2));
    if (hasROI){
      out_recowaveforms2->emplace_back(recob::OpWaveform(wf->TimeStamp(), wf->ChannelNumber(), rois));
      util::CreateAssn(*this, e, *out_recowaveforms2, wf, *recorawassn2);
    }
  }
  e.put(std::move(out_recowaveforms1),"wiener");
  e.put(std::move(out_recowaveforms2));
  e.put(std::move(recorawassn1),"wiener");
  e.put(std::move(recorawassn2));
}

DEFINE_ART_MODULE(pdsp::PDSPCalPhoton)
