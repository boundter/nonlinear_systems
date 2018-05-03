#ifndef __HELPER__
#define __HELPER__

#include <vector>

namespace nonlinear_systems {
/*! This helper class allows the extraction of the mean field for the calculation
 * of the average frequency, this way one observer can handle both general
 * systems and networks.
 */
template <typename mean_field_type>
class MeanFieldHelper {
  public:

    /*!
     *  Returns the phase of the given mean field. The phases are the elements,
     *  except for the first, so mean_field[1:]. It is overloaded to handle
     *  networks and flatten the phases.
     */
    std::vector<double> GetMeanFieldPhase(const mean_field_type& mean_field) {
      std::vector<double> phase;
      for (size_t i = 1; i < mean_field.size(); ++i) {
        phase.push_back(mean_field[i]);
      }
      return phase;
    }


    /*!
     *  Returns the order parameter of the given mean field. The order parameter
     *  is the first element in the mean field, so mean_field[0]. It is
     *  overloaded to handle networks and flatten the order parameters.
     */
    std::vector<double> GetOrderParameter(const mean_field_type& mean_field) {
      std::vector<double> order_parameter = {mean_field[0]};
      return order_parameter;
    }
};


// Template specialization for networks.
template<> std::vector<double> MeanFieldHelper<std::vector<std::vector<double> > >
::GetMeanFieldPhase(const std::vector<std::vector<double> >& mean_field) {
  std::vector<double> phase;
  for (unsigned int i = 0; i < mean_field.size(); ++i) {
    for (size_t j = 1; j < mean_field[i].size(); ++j) {
      phase.push_back(mean_field[i][j]);
    }
  }
  return phase;
}


template<> std::vector<double> MeanFieldHelper<std::vector<std::vector<double> > >
::GetOrderParameter(const std::vector<std::vector<double> >& mean_field) {
  std::vector<double> order_parameter;
  for (size_t i = 0; i < mean_field.size(); ++i) {
    order_parameter.push_back(mean_field[i][0]);
  }
  return order_parameter;
}
} // nonlinear_systems

#endif
