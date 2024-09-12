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

#include "TMath.h"

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
    void compute_pedestal(vector<short> const& waveform, double threshold);
  };

  PDSMeanRMSLZ::PDSMeanRMSLZ(fhicl::ParameterSet const& ps)
    : raw_adc_thresh{ps.get<short>("raw_adc_thresh")}
    , rms_thresh{ps.get<double>("rms_thresh")}
    , n_iter{ps.get<short>("n_iter")}
  {}

  void PDSMeanRMSLZ::Analyze(vector<short> const& waveform){

    //static int count = 0;
    
    double med = compute_median(waveform);
    mean_ = med;
    //cout<<count<<" median = "<<med<<endl;
    //++count;
//    vector<double> data(waveform.size());
//    for (size_t i = 0; i<waveform.size(); ++i){
//      data[i] = waveform[i] - med;
//    }
    compute_pedestal(waveform, raw_adc_thresh);
    //cout<<"mean = "<<mean_<<" rms = "<<rms_<<endl;
//    for (size_t i = 0; i<data.size(); ++i){
//      data[i] = waveform[i] - mean_;
//    }
    
    for (int i = 0; i<n_iter; ++i){
      compute_pedestal(waveform, rms_thresh * rms_);
//      for (size_t j = 0; j<data.size(); ++j){
//        data[j] = data[j] - mean_;
//      }
      //cout<<"i = "<<i<<" mean = "<<mean_<<" rms = "<<rms_<<endl;
    }

    return;
  }

  // based on https://github.com/SBNSoftware/sbndcode/blob/ad265e6e705ee7b67fc74865f2571263f2d3b582/sbndcode/Decoders/TPC/SBNDTPCDecoder_module.cc#L224
  double PDSMeanRMSLZ::compute_median(vector<short> const& adcs){
    size_t asiz = adcs.size();
    int imed = 0;
    double median;
    if (!asiz){
      return 0;
    }
    else{
      imed = TMath::Median(asiz, adcs.data()) + 0.01;  // add an offset to make sure the floor gets the right integer
      median = imed;

      // add in a correction suggested by David Adams, May 6, 2019
    
      size_t s1 = 0;
      size_t sm = 0;
      for (size_t i = 0; i < asiz; ++i) {
        if (adcs.at(i) < imed) s1++;
        if (adcs.at(i) == imed) sm++;
      }
      if (sm > 0) {
        float mcorr = (-0.5 + (0.5*(double) asiz - (double) s1)/ ((double) sm) );
        median += mcorr;
      }
      
      return median;
    }
  }

  void PDSMeanRMSLZ::compute_pedestal(vector<short> const& waveform, double threshold){
    vector<short> adcs;
    vector<bool> ignore(waveform.size(), false);
    for (size_t i = 0; i<waveform.size(); ++i){
      double adc_subtract = waveform[i] - mean_;
      if (abs(adc_subtract) > threshold){
        for (int j = int(i) - 50; j<= int(i) + 100; ++j){
          if (j>=0 and j<int(waveform.size())){
            ignore[j] = true;
          }
        }
      }
    }

    double mean = 0;
    rms_ = 0;
    int n = 0;
    for (size_t i = 0; i<waveform.size(); ++i){
      double adc_subtract = waveform[i] - mean_;
      if (!ignore[i]){
        adcs.push_back(waveform[i]);
        mean += adc_subtract;
        rms_ += adc_subtract*adc_subtract;
        ++n;
      }
    }
    if (n){
      mean /=n;
      rms_ = sqrt(rms_/n - mean*mean);
    }
    mean_ = compute_median(adcs);
    if (!mean_){
      mean_ = compute_median(waveform);
      rms_ = TMath::RMS(waveform.size(), waveform.data());
    }
  }
}

DEFINE_ART_CLASS_TOOL(pds::PDSMeanRMSLZ)
