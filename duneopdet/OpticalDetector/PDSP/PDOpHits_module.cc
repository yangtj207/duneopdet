////////////////////////////////////////////////////////////////////////
// Class:       PDOpHits
// Plugin Type: analyzer (Unknown Unknown)
// File:        PDOpHits_module.cc
//
// Generated at Mon Sep  9 16:43:29 2024 by Tingjun Yang using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "lardataobj/RecoBase/Hit.h"

#include "TTree.h"
#include <vector>

using namespace std;

namespace pdsp {
  class PDOpHits;
}


class pdsp::PDOpHits : public art::EDAnalyzer {
public:
  explicit PDOpHits(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  PDOpHits(PDOpHits const&) = delete;
  PDOpHits(PDOpHits&&) = delete;
  PDOpHits& operator=(PDOpHits const&) = delete;
  PDOpHits& operator=(PDOpHits&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Selected optional functions.
  void beginJob() override;

private:

  art::InputTag fOpHitLabel;
  // Declare member data here.
  TTree *pdtree;
  int run; // run number
  int event; // event number
  vector<int> channel;
  vector<float> peaktime;
  vector<float> rms;
  vector<float> peak;
  vector<float> integral;
  vector<int> multiplicity;
  vector<int> localindex;
  vector<int> starttick;
  vector<int> endtick;
};


pdsp::PDOpHits::PDOpHits(fhicl::ParameterSet const& p)
  : EDAnalyzer{p} ,
    fOpHitLabel(p.get<art::InputTag>("OpHitLabel"))
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void pdsp::PDOpHits::analyze(art::Event const& e)
{
  // Implementation of required member function here.
  run = e.run();
  event = e.id().event();
  channel.clear();
  peaktime.clear();
  rms.clear();
  peak.clear();
  integral.clear();
  multiplicity.clear();
  localindex.clear();
  starttick.clear();
  endtick.clear();
  auto hitListHandle = e.getHandle< std::vector<recob::Hit> >(fOpHitLabel);
  for (size_t i = 0; i<(*hitListHandle).size(); ++i){
    auto &hit = (*hitListHandle)[i];
    channel.push_back(hit.Channel());
    peaktime.push_back(hit.PeakTime());
    rms.push_back(hit.RMS());
    peak.push_back(hit.PeakAmplitude());
    integral.push_back(hit.Integral());
    multiplicity.push_back(hit.Multiplicity());
    localindex.push_back(hit.LocalIndex());
    starttick.push_back(hit.StartTick());
    endtick.push_back(hit.EndTick());
    pdtree->Fill();
  }
}

void pdsp::PDOpHits::beginJob()
{
  // Implementation of optional member function here.
  art::ServiceHandle<art::TFileService> tfs;
  pdtree = tfs->make<TTree>("pdtree", "PD hits");
  pdtree->Branch("run", &run);
  pdtree->Branch("event", &event);
  pdtree->Branch("channel", &channel);
  pdtree->Branch("peaktime", &peaktime);
  pdtree->Branch("rms", &rms);
  pdtree->Branch("peak", &peak);
  pdtree->Branch("integral", &integral);
  pdtree->Branch("multiplicity", &multiplicity);
  pdtree->Branch("localindex", &localindex);
  pdtree->Branch("starttick", &starttick);
  pdtree->Branch("endtick", &endtick);

}

DEFINE_ART_MODULE(pdsp::PDOpHits)
