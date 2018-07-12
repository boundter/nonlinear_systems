#ifndef __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_SYSTEM__
#define __REDUCED_M_KURAMOTO_SAKAGUCHI_WATANABE_STROGATZ_SYSTEM__

#include <vector>
#include <nonlinear_systems/odes/reduced_m_kuramoto_sakaguchi_watanabe_strogatz_ode.hpp>
#include <nonlinear_systems/systems/generic_network.hpp>

typedef std::vector<double> state_type;
typedef std::vector<state_type> network_type;
typedef std::vector<unsigned int> node_size_type;

namespace nonlinear_systems {
class ReducedMKuramotoSakaguchiWatanabeStrogatzSystem
  : public GenericNetwork<ReducedMKuramotoSakaguchiWatanabeStrogatzODE, double> {
  public:

};
} // nonlinear_systems
#endif
