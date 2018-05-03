#ifndef __POSITION_OBSERVER__
#define __POSITION_OBSERVER__

#include <vector>

namespace nonlinear_systems{
/*!
 *  Observe s the position of the system and saves it together with the time.
 */
template<typename state_type = std::vector<double> >
struct PositionObserver {
  std::vector<state_type>& _x;
  std::vector<double>& _t;

  /*!
   *  @param position the vector where the position will be saved
   *  @param t the vector where the time will be saved
   */
  PositionObserver(std::vector<state_type>& position, std::vector<double>& t)
    : _x(position), _t(t) {;}

  void operator()(const state_type& x, double t) {
    _x.push_back(x);
    _t.push_back(t);
  }
};
} // nonlinear_systems

#endif
