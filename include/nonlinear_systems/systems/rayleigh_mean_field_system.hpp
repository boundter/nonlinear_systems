#include <vector>
#include <nonlinear_systems/odes/rayleigh_mean_field_ode.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>

typedef std::vector<double> state_type;

namespace nonlinear_systems {
  class RayleighMeanFieldSystem
    : public GenericSystem<RayleighMeanFieldODE,std::vector<double>, 
    boost::numeric::odeint::runge_kutta4<state_type>>  {
      public:
        RayleighMeanFieldSystem(unsigned int N, double nonlinearity, 
            double coupling, const char* coupling_coordinate = "y")
          :GenericSystem(N, 2){
            state_type frequency;
            for (unsigned int i = 0; i < N; ++i) {
              frequency.push_back(1);
            }
            _coupling_coordinate = coupling_coordinate;
            if(_coupling_coordinate == "x") {
              ode = new RayleighMeanFieldODEX(N, frequency, nonlinearity, 
                  coupling);
            }
            else if(_coupling_coordinate == "y") {
              ode = new RayleighMeanFieldODEY(N, frequency, nonlinearity, 
                  coupling);
            }
            else {
              throw std::invalid_argument("Did not understand coupling_coordinate. Allowed values are 'x'and 'y'.");
            }
          }

      protected:
        const char* _coupling_coordinate;
    };

} // nonlinear_systems
