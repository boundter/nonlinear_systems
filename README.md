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
For the installation change into the source directory and do
```
cmake . && make install
```

#### Testing

The tests are integrated into cmake:
```
cmake . && make test
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

typedef std::vector<double> state_type;

class HarmonicOscillator: {
  public:
    double omega = 1.;

    HarmonicOscillator(void* params) {};

    void operator()(const state_type& x, state_type& dx, const double t) {
      dx[0] = x[1];
      dx[1] = -omega*omega*x[0];
    }
}

GenericSystem<HarmonicOscillator> system = GenericSystem<HarmonicOscillator>(1, 2, NULL);

double dt = 0.01;
unsigned int n_steps = 10000;
system.Integrate(dt, n_steps);
state_type x = system.GetPosition();
double t = system.GetTime();
std::cout << "t = " << t << " x = {" << x[0] << ", " << x[1] << "}" << std::endl;

```
