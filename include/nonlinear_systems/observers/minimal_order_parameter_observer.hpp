#ifndef __MINIMAL_ORDER_PARAMETER_OBSERVER__
#define __MINIMAL_ORDER_PARAMETER_OBSERVER__

#include <vector>
#include <nonlinear_systems/misc/helper.hpp>

namespace nonlinear_systems {
template<typename system_type, typename state_type = std::vector<double>, 
  typename mean_field_type = state_type>
/*!
 *  \brief Observes the minimal value of the order parameter of a system
 */
class MinimalOrderParameterObserver {
  public:
    system_type& _system;
    std::vector<double>& _minimal_order_parameter;

    /*!
     *  @param system the observed system.
     *  @param minimal_order_parameter the minimal observed order parameter.
     */
    MinimalOrderParameterObserver(system_type& system, 
        std::vector<double>& minimal_order_parameter)
   : _system(system), _minimal_order_parameter(minimal_order_parameter) {
      for (size_t i = 0; i < _minimal_order_parameter.size(); ++i) {
        _minimal_order_parameter[i] = 1e40;
      }
    }


    void operator()(const state_type& x, double t) {
      mean_field_type mean_field = _system.CalculateMeanFieldSpherical();
      CompareOrderParameter(mean_field);
    }


  protected:
    MeanFieldHelper<mean_field_type> _helper;

    void CompareOrderParameter(const mean_field_type& mean_field) {
      std::vector<double> order_parameter = _helper.GetOrderParameter(mean_field);
      for (size_t i = 0; i < order_parameter.size(); ++i) {
        if (order_parameter[i] < _minimal_order_parameter[i]) {
          _minimal_order_parameter[i] = order_parameter[i];
        }
      }
    }
};
} //nonlinear_systems

#endif
