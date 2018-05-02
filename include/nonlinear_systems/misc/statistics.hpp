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


// TODO: chekc size?
template<typename state_type=std::vector<double> >
void UpdateVariance(const state_type& x, unsigned int& step_number,
    state_type& average, state_type& variance) {
  state_type previous_average = average;
  UpdateAverage(x, step_number, average);
  if (step_number == 2) {
    variance = state_type(x.size(), 0);
  }
  else {
    // TODO: increase numerical stability
    for (size_t i = 0; i < variance.size(); ++i) {
      variance[i] *= (static_cast<double>(step_number) - 3.); // TODO: is this correct?
      variance[i] += (x[i]-previous_average[i])*(x[i]-average[i]);
      variance[i] /= (static_cast<double>(step_number) - 2.);
    }
  }
}
}// nonlinear_systems

#endif
