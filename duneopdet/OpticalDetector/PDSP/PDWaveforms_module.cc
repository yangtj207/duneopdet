////////////////////////////////////////////////////////////////////////
// Class:       PDWaveforms
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDWaveforms_module.cc
//
// Generated at Sun Sep  8 15:23:55 2024 by Tingjun Yang using cetskelgen
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

#include <iostream>

using namespace std;

namespace pdsp {
  class PDWaveforms;
}

class pdsp::PDWaveforms : public art::EDAnalyzer {
public:
  explicit PDWaveforms(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDWaveforms(PDWaveforms const&) = delete;
  PDWaveforms(PDWaveforms&&) = delete;
  PDWaveforms& operator=(PDWaveforms const&) = delete;
  PDWaveforms& operator=(PDWaveforms&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  // Declare member data here.
  art::ServiceHandle<art::TFileService> tfs;
  TH1D *hwf;

};


pdsp::PDWaveforms::PDWaveforms(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::PDWaveforms::analyze(art::Event const& e)
{

  auto extListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:external");
  if (!extListHandle){
    cout<<"extListHandle invalid"<<endl;
    return;
  }

  auto intListHandle = e.getHandle< std::vector<raw::OpDetWaveform> >("ssprawdecoder:internal");
  if (!intListHandle){
    cout<<"intListHandle invalid"<<endl;
    return;
  }

  art::FindOneP<recob::OpWaveform> extip(extListHandle, e, "calpdext");
  art::FindOneP<recob::OpWaveform> extwn(extListHandle, e, "calpdext:wiener");
  
  art::FindOneP<recob::OpWaveform> intip(intListHandle, e, "calpdint");
  art::FindOneP<recob::OpWaveform> intwn(intListHandle, e, "calpdint:wiener");

  for (size_t i = 0; i<extListHandle->size(); ++i){
    auto &wf = (*extListHandle)[i];
    int daqch = wf.ChannelNumber();
    hwf = tfs->make<TH1D>(Form("ext_%d",daqch), Form("ext, ch %d", daqch), 2000, 0, 2000);
    for (unsigned short j = 0; j<wf.Waveform().size(); ++j){
      hwf->SetBinContent(j+1, wf.Waveform()[j]);
    }
    if (extwn.at(i).isAvailable()){
      auto &op = extwn.at(i);
      hwf = tfs->make<TH1D>(Form("extwn_%d",daqch), Form("extwn, ch %d", daqch), 2000, 0, 2000);
      for (unsigned short j = 0; j<op->Signal().size(); ++j){
        hwf->SetBinContent(j+1, op->Signal()[j]);
      }
    }
    if (extip.at(i).isAvailable()){
      auto &op = extip.at(i);
      hwf = tfs->make<TH1D>(Form("extip_%d",daqch), Form("extip, ch %d", daqch), 2000, 0, 2000);
      for (unsigned short j = 0; j<op->Signal().size(); ++j){
        hwf->SetBinContent(j+1, op->Signal()[j]);
      }
    }
  }

  for (size_t i = 0; i<intListHandle->size(); ++i){
    auto &wf = (*intListHandle)[i];
    int daqch = wf.ChannelNumber();
    hwf = tfs->make<TH1D>(Form("int_%d_%d",daqch,int(i)), Form("int, ch %d", daqch), 2000, 0, 2000);
    for (unsigned short j = 0; j<wf.Waveform().size(); ++j){
      hwf->SetBinContent(j+1, wf.Waveform()[j]);
    }

    if (intwn.at(i).isAvailable()){
      auto &op = intwn.at(i);
      hwf = tfs->make<TH1D>(Form("intwn_%d_%d",daqch,int(i)), Form("intwn, ch %d", daqch), 2000, 0, 2000);
      for (unsigned short j = 0; j<op->Signal().size(); ++j){
        hwf->SetBinContent(j+1, op->Signal()[j]);
      }
    }

    if (intip.at(i).isAvailable()){
      auto &op = intip.at(i);
      hwf = tfs->make<TH1D>(Form("intip_%d_%d",daqch,int(i)), Form("intip, ch %d", daqch), 2000, 0, 2000);
      for (unsigned short j = 0; j<op->Signal().size(); ++j){
        hwf->SetBinContent(j+1, op->Signal()[j]);
      }
    }
  }

}

void pdsp::PDWaveforms::beginJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(pdsp::PDWaveforms)
