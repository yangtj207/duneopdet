////////////////////////////////////////////////////////////////////////
// Class:       PDSMeanRMSLZ
// Plugin Type: tool
// File:        PDSMeanRMSLZ_tool.cc
// Author:      tjyang@fnal.gov
// 
// Calculate pedestal and RMS of a photon detector waveform
// Based on code implemented by Laura Zamlli
// Ref: https://github.com/dune-lardon/lardon/blob/main/pedestals.py
// 
////////////////////////////////////////////////////////////////////////

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "PDSMeanRMS.h"

#include "TH1D.h"

using namespace std;

namespace pds{

  class PDSMeanRMSLZ : public PDSMeanRMS {
  public:
    explicit PDSMeanRMSLZ(fhicl::ParameterSet const& ps);

    void Analyze(vector<short> const& waveform) override;
    double GetMean() override {return mean_;}
    double GetRMS() override {return rms_;}

  private:
    double mean_;
    double rms_;
    short raw_adc_thresh;
    double rms_thresh;
    short n_iter;
    double compute_median(vector<short> const& scores);
    void compute_pedestal(vector<double> const& waveform, double threshold);
  };

  PDSMeanRMSLZ::PDSMeanRMSLZ(fhicl::ParameterSet const& ps)
    : raw_adc_thresh{ps.get<short>("raw_adc_thresh")}
    , rms_thresh{ps.get<double>("rms_thresh")}
    , n_iter{ps.get<short>("n_iter")}
  {}

  void PDSMeanRMSLZ::Analyze(vector<short> const& waveform){

    double med = compute_median(waveform);
    vector<double> data(waveform.size());
    for (size_t i = 0; i<waveform.size(); ++i){
      data[i] = waveform[i] - med;
    }
    compute_pedestal(data, raw_adc_thresh);
    for (size_t i = 0; i<data.size(); ++i){
      data[i] = data[i] - mean_;
    }
    
    for (int i = 0; i<n_iter; ++i){
      compute_pedestal(data, rms_thresh * rms_);
      for (size_t j = 0; j<data.size(); ++j){
        data[j] = data[j] - mean_;
      }
    }

    return;
  }

  double PDSMeanRMSLZ::compute_median(vector<short> const& scores){
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
  
  void PDSMeanRMSLZ::compute_pedestal(vector<double> const& waveform, double threshold){
    mean_ = 0;
    rms_ = 0;
    int n = 0;
    for (auto const & adc : waveform){
      if (adc < threshold){
        mean_ += adc;
        rms_ += adc*adc;
        ++n;
      }
    }
    if (n){
      mean_ /=n;
      rms_ = sqrt(rms_/n - mean*mean);
    }
  }
}

DEFINE_ART_CLASS_TOOL(pds::PDSMeanRMSLZ)
