#include <vector>

namespace pds{
  
  class PDSMeanRMS {
  public:

    virtual ~PDSMeanRMS() noexcept = default;
    virtual void Analyze(std::vector<short> const& waveform) = 0;
    virtual double GetMean() = 0;
    virtual double GetRMS() = 0;
  };
}
