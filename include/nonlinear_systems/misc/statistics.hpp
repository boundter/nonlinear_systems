#ifndef __STATISTICS__
#define __STATISTICS__

#include <vector>


namespace nonlinear_systems{
/*!
 *  Sample a given distribution.
 */
template<typename state_type=std::vector<double>,typename preci=double>
state_type SampleDistribution(size_t number_samples, 
    std::function<preci()>* distribution) {
  state_type samples;
  for (size_t i = 0; i < number_samples; ++i) {
    samples.push_back((*distribution)());
  }
  return samples;
}


// TODO: check size?
template<typename state_type=std::vector<double> > 
void UpdateAverage(const state_type& x, unsigned int& step_number, 
    state_type& average) {
  for (size_t i = 0; i < average.size(); ++i) {
    average[i] += (x[i] - average[i])/static_cast<double>(step_number);
  } 
  step_number += 1;
}
}// nonlinear_systems

#endif
