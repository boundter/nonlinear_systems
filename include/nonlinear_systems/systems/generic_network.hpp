#ifndef __GENERIC_NETWORK__
#define __GENERIC_NETWORK__

#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>

#include <iostream>

// TODO: where should be virtual?

namespace nonlinear_systems {
template<typename ode_type, typename precision_type = double, 
  typename stepper_type = 
    boost::numeric::odeint::runge_kutta4<std::vector<precision_type> > >
class GenericNetwork: protected GenericSystem<ode_type, 
  std::vector<precision_type>, stepper_type> {
  public:
  typedef std::vector<precision_type> state_type;
  typedef std::vector<unsigned int> node_size_type;
  typedef std::vector<state_type> matrix_type;
    
    GenericNetwork(node_size_type node_sizes, unsigned int dimension,
        void* parameters)
      : GenericSystem<ode_type, state_type, stepper_type>(0, dimension, 
          parameters) {
        _node_indices = CalculateNodeIndices(node_sizes);
        this->Resize(_node_indices.back());
    }


    void SetState(const state_type& new_state) {
      this->SetPosition(new_state);
    }
 

    state_type GetState() {
      return GenericSystem<ode_type, state_type, stepper_type>::GetPosition();
    }


    matrix_type GetNodes() {
      matrix_type nodes;
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        nodes.push_back(state_type());
        for (unsigned int j = _node_indices[i]; j < _node_indices[i+1]; ++j) {
          nodes.back().push_back(this->_x[j]);
        }
      }
      return nodes;
    }

  
    node_size_type GetNodeIndices() {
      return _node_indices;
    }


    double GetTime() {
      GenericSystem<ode_type, state_type, stepper_type>::GetTime();  
    }


    double SetParameters(void* parameters) {
      GenericSystem<ode_type, state_type, stepper_type>::SetParameters(parameters);
    }


    template<typename observer_type = boost::numeric::odeint::null_observer>
    void Integrate(double dt, unsigned int number_steps, 
        observer_type observer = observer_type()) {
      GenericSystem<ode_type, state_type, stepper_type>::template
        Integrate<observer_type>(dt, number_steps, observer);
    }


    matrix_type CalculateMeanField() {
      matrix_type mean_field;
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        mean_field.push_back(state_type(this->_d));
        unsigned int number_oscillators = _node_indices[i+1] - _node_indices[i];
        for (unsigned int k = 0; k < this->_d; ++k) {
          for (unsigned int j = this->_d*_node_indices[i] + k; 
              j < this->_d*_node_indices[i+1]; j += this->_d) {
            mean_field.back()[k] += this->_x[j];
          }
          mean_field.back()[k] /= static_cast<double>(number_oscillators);
        }
      }
      return mean_field;
    }

  protected:
    node_size_type _node_indices;
    
    node_size_type CalculateNodeIndices(const node_size_type& node_sizes) {
      node_size_type node_indices = {0};
      unsigned int offset = 0;;
      for (size_t i = 0; i < node_sizes.size(); ++i) {
        node_indices.push_back(offset + node_sizes[i]);
        offset += node_sizes[i];
      }
      return node_indices;
    }
};
} // nonlinear_systems

#endif
