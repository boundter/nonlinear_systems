#ifndef __STATISTICS__
#define __STATISTICS__

#include <vector>


namespace nonlinear_systems{

template<typename state_type=std::vector<double>,typename preci=double>
state_type SampleDistribution(size_t number_samples, 
    std::function<preci()>* distribution) {
  state_type samples;
  for (size_t i = 0; i < number_samples; ++i) {
    samples.push_back((*distribution)());
  }
  return samples;
}


}// nonlinear_systems

#endif
