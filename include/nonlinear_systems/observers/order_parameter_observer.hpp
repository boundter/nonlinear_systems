#ifndef __ORDER_PARAMETER_OBSERVER__
#define __ORDER_PARAMETER_OBSERVER__

#include <vector>
#include <nonlinear_systems/observers/statistics_observer.hpp>
#include <nonlinear_systems/misc/helper.hpp>

namespace nonlinear_systems {
template<typename system_type, typename state_type = std::vector<double>,
  typename mean_field_type = state_type >
class VarianceOrderParameterObserver {
  public:
    system_type& _system;

    VarianceOrderParameterObserver(system_type& system, 
        state_type& average_order_parameter, state_type& variance_order_parameter)
      : _system(system) {
        _variance_observer = 
          std::shared_ptr<VarianceObserver<std::vector<double>>>(
              new VarianceObserver<std::vector<double>>(average_order_parameter,
                variance_order_parameter));
      }


      void operator()(const state_type& x, double t) {
        mean_field_type mean_field = _system.CalculateMeanFieldSpherical();
        std::vector<double> order_parameter = GetOrderParameter(mean_field);
        _variance_observer->operator()(order_parameter, t);
      }

  protected:
    std::shared_ptr<VarianceObserver<std::vector<double>>> _variance_observer;
    MeanFieldHelper<mean_field_type> _mean_field_helper;

    std::vector<double> GetOrderParameter(const mean_field_type& mean_field) {
      return _mean_field_helper.GetOrderParameter(mean_field);
    }
};
} // nonlinear_systems
#endif
