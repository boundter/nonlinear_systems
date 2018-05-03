#ifndef __STATISTICS_OBSERVER__
#define __STATISTICS_OBSERVER__

#include <vector>
#include <nonlinear_systems/misc/statistics.hpp>

namespace nonlinear_systems {
/*!
 * Observes the moving average of the given quantity.
 */
template<typename state_type = std::vector<double> >
class AverageObserver {
  public:
    std::vector<double>& _average;
    unsigned int _step_number;

    /*!
     *  @param average the vector in wich to save the moving average.
     */
    AverageObserver(std::vector<double>& average)
    : _average(average) {
      _step_number = 1;
    }


    void operator()(const state_type& x, double t) {
      UpdateAverage(x, _step_number, _average);
    }
};


/*!
 *  Observes the moving average and variance of the given quantity. The variance
 *  is define as 
 *  \f[ Var(x) = \frac{\sum_{j=1}^{N} (x_j - Avg(x))}{N-1} \f]
 */
template<typename state_type = std::vector<double> >
class VarianceObserver {
  public:
    std::vector<double>& _average;
    std::vector<double>& _variance;
    unsigned int _step_number;

    VarianceObserver(std::vector<double>& average, std::vector<double>& variance)
    :_average(average), _variance(variance) {
      _step_number = 1;
    }

    
    void operator()(const state_type& x, double t) {
      UpdateSumOfSquares(x, _step_number, _average, _sum_of_squares);
      if (_step_number > 2) {
        for (size_t i = 0; i < _sum_of_squares.size(); ++i) {
          // divide by _step_number - 2, instead of _step_number - 1, because it
          // was already increased in UpdateSumOfSquares
          _variance[i] = _sum_of_squares[i]/static_cast<double>(_step_number - 2);
        }
      }
    }

  private:
    std::vector<double> _sum_of_squares;
};
} // nonlinear_systems

#endif
