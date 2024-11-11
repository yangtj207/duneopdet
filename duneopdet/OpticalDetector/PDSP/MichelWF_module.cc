////////////////////////////////////////////////////////////////////////
// Class:       MichelWF
// Plugin Type: analyzer (Unknown Unknown)
// File:        MichelWF_module.cc
//
// Generated at Sat Sep 14 19:07:15 2024 by Tingjun Yang using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardataobj/RecoBase/OpWaveform.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF1Convolution.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;

double fitf(double *x, double *par){

  double fitval = par[0]/(3*par[1]/20)*exp(-x[0]/(3*par[1]/20))+(1-par[0])/(3*par[2]/20)*exp(-x[0]/(3*par[2]/20));
                                    
  return fitval;
}

double fitf2(double *x, double *par){

  double fitval = par[0]*exp(-(x[0]-par[1])*(x[0]-par[1])/par[2]/par[2]);
                                    
  return fitval;
}

namespace pdsd {
  class MichelWF;
}

class pdsd::MichelWF : public art::EDAnalyzer {
public:
  explicit MichelWF(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MichelWF(MichelWF const&) = delete;
  MichelWF(MichelWF&&) = delete;
  MichelWF& operator=(MichelWF const&) = delete;
  MichelWF& operator=(MichelWF&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;
  void endJob() override;

private:

  // Declare member data here.
  TH1D *hwf;
  TH1D *hdt;
  TH2D *hdq;
  TF1 *fun;
  vector<double> pred;
  TH1D *hwf1;
  TH1D *hwf2;
  TCanvas *can;
};


pdsd::MichelWF::MichelWF(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
  TF1 *exps = new TF1("exps", fitf, 0, 2000, 3);
  TF1 *fgaus = new TF1("fgaus", fitf2, 0, 2000, 3);
  TF1Convolution *f_conv = new TF1Convolution(exps,fgaus,0,2000,true);
  f_conv->SetRange(0.,2000.);
  f_conv->SetNofPointsFFT(10000);
  //cout<<f_conv->GetNpar()<<endl;
  TF1   *fun = new TF1("fit",*f_conv, 0., 2000., f_conv->GetNpar());
  fun->SetParameters(0.09179, 0.3174, 1317.44, 2.498e4/11681.3, 66.94, 13.14);
  for (int i = 0; i<2000; ++i) pred.push_back(fun->Eval(i+0.5));
  hwf1 = new TH1D("hwf1","hwf1",2000,0,2000);
  hwf2 = new TH1D("hwf2","hwf2",2000,0,2000);
}

void pdsd::MichelWF::analyze(art::Event const& e)
{
  // Implementation of required member function here.

//  auto extListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:external");
//  if (!extListHandle){
//    cout<<"extListHandle invalid"<<endl;
//    return;
//  }

  auto intListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:internal");
  if (!intListHandle){
    cout<<"intListHandle invalid"<<endl;
    return;
  }

  art::FindOneP<recob::OpWaveform> intip(intListHandle, e, "calpdint");
  art::FindOneP<recob::OpWaveform> intwn(intListHandle, e, "calpdint:wiener");

  //auto &wfext = (*extListHandle)[0];
  vector<bool> used(intListHandle->size(),false);
  vector<double> sumwf(2000, 0);
  for (size_t i = 0; i<intListHandle->size(); ++i){
    auto &wf = (*intListHandle)[i];
    int daqch = wf.ChannelNumber();
    if ((daqch>=132&&daqch<=143)||
        (daqch>=264&&daqch<=275)){
      if (intip.at(i).isAvailable()){
        auto &op = intip.at(i);
//        int maxbin = -1;
//        double maxpe = -1;
//        for (unsigned short j = 0; j<op->Signal().size(); ++j){
//          if (op->Signal()[j]>maxpe){
//            maxbin = j;
//            maxpe = op->Signal()[j];
//          }
//        }
        for (unsigned short j = 0; j<op->Signal().size(); ++j){
//          int bin = j - maxbin + 70;
//          if (bin>0 && bin<2000 && maxbin < 100){
          hwf->SetBinContent(j+1, hwf->GetBinContent(j+1)+op->Signal()[j]);
            //          }
        }
      }
      if (!used[i]){
        if (intip.at(i).isAvailable()){
          auto &op = intip.at(i);
          for (unsigned short j = 0; j<op->Signal().size(); ++j){
            sumwf[j] = op->Signal()[j];
          }
        }
        used[i] = true;
      }
      else{
        continue;
      }
      for (size_t j = 0; j<intListHandle->size(); ++j){
        auto &wf1 = (*intListHandle)[j];
        int daqch1 = wf1.ChannelNumber();
        if (j==i) continue;
        if ((daqch1>=132&&daqch1<=143)||
            (daqch1>=264&&daqch1<=275)){
          //cout<<wf.TimeStamp() - wf1.TimeStamp()<<endl;
          hdt->Fill(wf.TimeStamp() - wf1.TimeStamp());
          if (abs(wf.TimeStamp() - wf1.TimeStamp())<0.2){
            used[j] = true;
            if (intip.at(j).isAvailable()){
              auto &op = intip.at(j);
              for (unsigned short k = 0; k<op->Signal().size(); ++k){
                sumwf[k] += op->Signal()[k];
              }
            }
          }
        }
      }
      double peak = 0;
      for (size_t j = 0; j<100; ++j){
        if (sumwf[j] > peak){
          peak = sumwf[j];
        }
      }
      //cout<<peak<<endl;
      for (size_t j = 30; j<2000; ++j){
        //cout<<sumwf[j]<<" "<<peak*pred[j]<<endl;
        hdq->Fill(j, sumwf[j] - peak*pred[j]);
      }
      for (int j = 0; j<2000; ++j){
        hwf1->SetBinContent(j+1, sumwf[j]);
        hwf2->SetBinContent(j+1, peak*pred[j]);
      }
      hwf1->SetTitle(Form("Run %d Event %d %d", e.run(), e.id().event(), int(i)));
      hwf1->Draw();
      hwf2->SetLineColor(2);
      hwf2->Draw("same");
      can->Print("wfs.ps");
      for (size_t j = 0; j<2000; ++j) sumwf[j] = 0;
    }
  }

}

void pdsd::MichelWF::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  hwf = tfs->make<TH1D>("intwf", "waveform", 2000, 0, 2000);
  hdt = tfs->make<TH1D>("hdt", "delta t", 1000, -10, 10);
  hdq = tfs->make<TH2D>("hdq", "delta q", 2000, 0, 2000, 100, -10,10);
  can = new TCanvas("can", "can");
  can->Print("wfs.ps[");
  
}

void pdsd::MichelWF::endJob()
{

  can->Print("wfs.ps]");
  
}

DEFINE_ART_MODULE(pdsd::MichelWF)
