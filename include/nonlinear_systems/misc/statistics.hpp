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


/*!
 *  Updates the moving average.
 *
 *  @param x the current state.
 *  @param step_number the number of previous steps/ data points. It will be
 *  incremented automatically.
 *  @param average the average of the previous state. It will be overwritten
 *  with the new average.
 */
template<typename state_type=std::vector<double> > 
void UpdateAverage(const state_type& x, unsigned int& step_number, 
    state_type& average) {
  for (size_t i = 0; i < average.size(); ++i) {
    average[i] += (x[i] - average[i])/static_cast<double>(step_number);
  } 
  step_number += 1;
}


// TODO: chekc size?
/*!
 *  Updates the moving average and sum of squares, needed for the variance. The
 *  variance can be calculated by dividing the sum of squares by the number of
 *  datapoints - 1, but the sum of squares is numerically more stable than the
 *  variance.
 *
 *  @param x the current state.
 *  @param step_number the number of tspes/ data points. It will be incremented
 *  automatically.
 *  @param average the average of the previous state. It will be overwritten
 *  with the new average.
 *  @param sum_of_squares the sum of suqares of the previous state. It will be
 *  overwritten with the new average.
 */
template<typename state_type=std::vector<double> >
void UpdateSumOfSquares(const state_type& x, unsigned int& step_number,
    state_type& average, state_type& sum_of_squares) {
  state_type previous_average = average;
  UpdateAverage(x, step_number, average);
  // for the first step there is no sum of squares, this corresponds to
  // step_number == 2, because it was already incremented in UpdateAverage
  if (step_number == 2) {
    sum_of_squares = state_type(x.size(), 0);
  }
  else {
    for (size_t i = 0; i < sum_of_squares.size(); ++i) {
      sum_of_squares[i] += (x[i]-previous_average[i])*(x[i]-average[i]);
    }
  }
}
}// nonlinear_systems

#endif
