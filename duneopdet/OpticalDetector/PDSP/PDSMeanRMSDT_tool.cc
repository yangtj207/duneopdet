////////////////////////////////////////////////////////////////////////
// Class:       PDSMeanRMSDT
// Plugin Type: tool
// File:        PDSMeanRMSDT_tool.cc
// Author:      tjyang@fnal.gov
// 
// Calculate pedestal and RMS of a photon detector waveform
// Based on code implemented by Dante Totani
// 
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "PDSMeanRMS.h"

#include "TH1D.h"

namespace pds{

  class PDSMeanRMSDT : public PDSMeanRMS {
  public:
    explicit PDSMeanRMSDT(fhicl::ParameterSet const& ps);

    void Analyze(std::vector<short> const& waveform) override;
    double GetMean() override {return mean_;}
    double GetRMS() override {return rms_;}

  private:
    double mean_;
    double rms_;
  };

  PDSMeanRMSDT::PDSMeanRMSDT(fhicl::ParameterSet const& ps)
  {}

  void PDSMeanRMSDT::Analyze(std::vector<short> const& waveform){

    TH1D *basehelp= new TH1D("basehelp","basehelp",2000, 1300,1800);
    for(size_t j=0; j<waveform.size(); j++){
      if(j<1000){
        basehelp->Fill(waveform[j]);
      }
    }
    int basebinmax = basehelp->GetMaximumBin();
    mean_ = basehelp->GetXaxis()->GetBinCenter(basebinmax);
    rms_  = basehelp->GetRMS();
    basehelp->Delete();   
    return;
  }
}

DEFINE_ART_CLASS_TOOL(pds::PDSMeanRMSDT)
