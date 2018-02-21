# Nonlinear Systems

Nonlinear Systems contains a generic interface for integrating ordinary differential equations. This reduces the repetitive work of integrating such systems.

## Getting Started
This is (for now) a header-only library, so there is no need to install anything.

### Prerequisites
- boost (at least version 1.61, this is the oldest version it was tested with)
- for installation: cmake (lowest tested version is 3.0)
- for documentations: doxygen

### Installing (optional)

Clone the repository:
```
git clone https://github.com/boundter/nonlinear_systems
```
For the installation change into the source directory and do:
```
cmake . && make install
```

To use the library in a program it has to be gnerated with cmake in the ```CMakeLists.txt``` it can be found like any other package:
```
find_package(nonliear_systems REQUIRED)
```
It can than be linked to any target:
```
target_link_libraries(target nonlinear_systems)
```
When including it like this the ```Boost_INCLUDE_DIRS``` has to be set and included. An example CMakeLists could look like this
```
cmake_minimum_required(VERSION 3.0)
project(TEST_APP CXX)

find_package(nonlinear_systems REQUIRED)
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(test_app test_app.cpp)
target_link_libraries(test_app nonlinear_systems)
```


#### Testing

The tests are integrated into cmake:
```
cmake . && make && make test
```

### Documentation

The documentation can be generated with doxygen.
```
cd docs && doxygen
```

The html-doc can then be found in ``` html/index.html ```.

## Example
Say you want to integrate the harmonic oscillator and get the position after some time. Assuming the files were installed it will be
```
#include <vector>
#include <iostream>
#include <nonlinear_systems/systems/generic_system.hpp>

using namespace std;
using namespace nonlinear_systems;

typedef vector<double> state_type;

// Define the ODE to solve with frequency omega = 1
class HarmonicOscillator {
  public:
    double omega = 1.;

    HarmonicOscillator(void* params) {};

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega*omega*x[0];
    }
};

int main() {
  GenericSystem<HarmonicOscillator> system = GenericSystem<HarmonicOscillator>(1, 2, NULL);

  // Set the initial condition
  state_type x_0 {1., 0.};
  system.SetPosition(x_0);

  // Integrate the system with timestep 0.01 for 10000 steps
  double dt = 0.01;
  unsigned int n_steps = 10000;
  system.Integrate(dt, n_steps);

  // Get the position
  state_type x = system.GetPosition();
  double t = system.GetTime();
  cout << "t = " << t << ", x = {" << x[0] << ", " << x[1] << "}" << endl;
}
```
