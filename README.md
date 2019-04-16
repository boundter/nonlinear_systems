# Nonlinear Systems

Nonlinear Systems contains a generic interface for integrating ordinary differential equations (ODEs). This reduces the repetitive work of integrating such systems.

## Used tools
- C++
- CMake to build, install and test
- Doxygen for documentation
- Boost
  - odeint for integration of the ODEs
  - program_options for adding command line arguments
  - test for executing the unit-tests
- [NLopt](https://github.com/stevengj/nlopt) (may become voluntary in the future)

## Whats In There?

nonlinear_systems is a framework to increase the readability of numerical integration of ODEs. A nice side effect is the reduction of bugs. The most important part is the ```GenericSystem``` in ```<nonlinear_systems/systems/generic_system.hpp>```. This class acts as a wrapper around the ODE and keeps track of the current state of the system. It takes the ODE as a type, where the ```operator()``` of the ODE calculates the ODE for given ```x``` and ```t```. The interaction with the system is then handled by ```GenericSystem```.

There are simple methods to get the current time or state, as well as more complicated ones to calculate the period for a given poincare surface. For a system of multiple units, there is also the posibility to calculate the mean field. The observation during the integration are handled by [observers for odeint](https://www.boost.org/doc/libs/1_66_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/short_example.html). Some generic ones are included in ```nonlinear_systems/observers```. For more complicated systems, where different methods are needed, there is the possibility to inherit from ```GenericSystem```.

A similar framework for networks of ODEs is also available in ```GenericNetwork```. The functions are similar and allow a conversion between thinking in terms of 2d-vectors, as well as a flattened representation.

Finally there are some ODEs prepacked, like the Kuramoto-model or the harmonic oscillator. Everything mentioned is unit tested, although the coverage is not quite 100%.

## Getting Started
This is a header-only library, so there is no need to install anything. It is a wrapper around the odeint-boost package

### Prerequisites
- boost (at least version 1.61, this is the oldest version it was tested with)
- for installation: cmake (lowest tested version is 3.0)
- for documentations: doxygen
- NLopt (at least version 2.5)

### Installing (optional)

Clone the repository:
```
git clone https://github.com/boundter/nonlinear_systems
```
For the installation change into the source directory and do:
```
cmake . && make install
```

To use the library in a program it has to be build with cmake. It can then be found like any other package.
```
find_package(nonliear_systems REQUIRED)
```
To link the library just use
```
target_link_libraries(target nonlinear_systems)
```
When including it like this, the ```Boost_INCLUDE_DIRS``` and ```NLOPT_INCLUDE_DIRS``` has to be set and included. An example CMakeLists could look like this
```
cmake_minimum_required(VERSION 3.0)
project(TEST_APP CXX)

find_package(nonlinear_systems REQUIRED)
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

find_package(NLopt 2.5.0)
include_directories(${NLOPT_INCLUDE_DIRS})

add_executable(test_app test_app.cpp)
target_link_libraries(test_app nonlinear_systems)
```


#### Testing

The tests are integrated into cmake and can be build by setting the variable `NL_BUILD_TESTS` to `ON`:
```
cmake -DNL_BUILD_TESTS=ON . && make && make test
```

### Documentation

The documentation can be generated with doxygen.
```
cd docs && doxygen
```

The html-doc can then be found in ``` html/index.html ```.

## Example
Say you want to integrate the harmonic oscillator and calculate the position.
```
#include <iostream>
#include <vector>
#include <nonlinear_systems/systems/generic_system.hpp>
#include <nonlinear_systems/observers/position_observer.hpp>

namespace nl = nonlinear_systems;

typedef  std::vector<double> state_type;

class HarmonicOscillatorODE {
  public:
    double _omega;

    HarmonicOscillatorODE(void* params) {
      _omega = reinterpret_cast<double*>(params)[0];
    }

    void operator()(const state_type& x, state_type& dx, double t) {
      dx[0] = x[1];
      dx[1] = -_omega*_omega*x[0];
    }
};

int main() {
  double omega[] = {1.};  // define the frequency
  auto HarmonicOsc = nl::GenericSystem<HarmonicOscillatorODE>(1, 2, omega);
  state_type initial = {0., 1.};
  HarmonicOsc.SetPosition(initial);  // set the initial condition
  std::vector<state_type> x; std::vector<double> t;  // container for the pos.
  double dt = 0.1; int n_steps = 10;
  HarmonicOsc.Integrate(dt, n_steps, nl::PositionObserver<state_type>(x, t));
  for (size_t i = 0; i < x.size(); ++i) {
    std::cout << t[i] << ": (" <<
      x[i][0] << ", " << x[i][1] << ")" << std::endl;
  }
}
```
The output of this program is
```
0: (0, 1)
0.1: (0.0998333, 0.995004)
0.2: (0.198669, 0.980067)
0.3: (0.29552, 0.955337)
0.4: (0.389418, 0.921061)
0.5: (0.479425, 0.877583)
0.6: (0.564642, 0.825336)
0.7: (0.644217, 0.764843)
0.8: (0.717356, 0.696707)
0.9: (0.783326, 0.621611)
1: (0.84147, 0.540303)
```
