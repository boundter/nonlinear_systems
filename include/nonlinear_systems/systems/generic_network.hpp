#ifndef __GENERIC_NETWORK__
#define __GENERIC_NETWORK__

#include <vector>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <nonlinear_systems/systems/generic_system.hpp>

#include <iostream>

// TODO: where should be virtual?

namespace nonlinear_systems {
/*!
 *  This class is a wrapper for the GenericSystem for the use with a network.
 *  A network consists of linked nodes, where each node can consists of
 *  multiple oscillators. There is no limitation as to the type of oscillator,
 *  it can be a phase oscillator or a general limit cycle one. For now it is
 *  limited to a use of oscillators of the same dimensionality.
 */
template<typename ode_type, typename precision_type = double,
  typename stepper_type =
    boost::numeric::odeint::runge_kutta4<std::vector<precision_type> > >
class GenericNetwork: protected GenericSystem<ode_type,
  std::vector<precision_type>, stepper_type> {
  public:
  typedef std::vector<precision_type> state_type;
  typedef std::vector<unsigned int> node_size_type;
  typedef std::vector<state_type> matrix_type;

    /*!
     * The network is initialized to a zero state.
     *
     * @param node_sizes a vector containing the size of every single node
     * @param dimension the dimensionality of the oscillators
     * @param parameters pointer to parameters for the ODE
     */
    GenericNetwork(node_size_type node_sizes, unsigned int dimension,
        void* parameters)
      : GenericSystem<ode_type, state_type, stepper_type>(0, dimension,
          parameters) {
        _node_indices = CalculateNodeIndices(node_sizes);
        _node_sizes = node_sizes;
        this->Resize(CalculateNumberOscillators(_node_sizes));
    }


    /*!
     * Sets the state of the system using a flattened representation of the form
     * state = {node_1x_1, node_1x_2, ..., node_2x_1, ....}.
     */
    void SetPosition(const state_type& new_state) {
      GenericSystem<ode_type, state_type, stepper_type>::SetPosition(new_state);
    }


    /*!
     * Gets the state in a flattened representation of the from
     * state = {node_1x_1, node_1x_2, ..., node_2x_1, ....}.
     */
    state_type GetPosition() {
      return GenericSystem<ode_type, state_type, stepper_type>::GetPosition();
    }


    /*!
     *  Gets the derivative in a flattened representation.
     */
    state_type GetDerivative() {
      return GenericSystem<ode_type, state_type, stepper_type>::GetDerivative();
    }


    /*!
     * \brief Return the position in the state space in phases for all elements.
     *
     * Calculates the coordinates on a sphere of the same dimension as
     * the phase space. If the dimension is 1, the corrdinates will be wrapped
     * around the unit circle as phases, otherise the first coordinate of every
     * element is the radius and the later ones are the phases.
     * Careful: in 3-d this is not the same as spherical coordinates with polar
     * angle and azimuth!
     */
    state_type GetPositionSpherical() {
      return GenericSystem<ode_type, state_type, stepper_type>::GetPositionSpherical();
    }


    /*!
     *  Gets the state as avector of vector representation
     *  state = {{node_1x_1, node_1x_2, ....}, {node_2x_1, ...}, ...}.
     */
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


    /*!
     * Gets the indices of the beginning of every new node + (the last index + 1)
     * of the flattened representation.
     */
    node_size_type GetNodeIndices() {
      return _node_indices;
    }


    /*!
     * Gets the time of the system.
     */
    double GetTime() {
      GenericSystem<ode_type, state_type, stepper_type>::GetTime();
    }


    /*!
     * Sets the parameters of the ODE.
     */
    void SetParameters(void* parameters) {
      GenericSystem<ode_type, state_type, stepper_type>::SetParameters(parameters);
    }


    /*!
     * Integrate the system in time.
     *
     * @param dt timestep
     * @param number_steps number of steps in time
     * @param observer observer during the integration
     */
    template<typename observer_type = boost::numeric::odeint::null_observer>
    void Integrate(double dt, unsigned int number_steps,
        observer_type observer = observer_type()) {
      GenericSystem<ode_type, state_type, stepper_type>::template
        Integrate<observer_type>(dt, number_steps, observer);
    }


    /*!
     * Calculate the mean field as the average of the position of the
     * oscillators in every node.
     */
    virtual matrix_type CalculateMeanField() {
      matrix_type mean_field;
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        mean_field.push_back(GenericSystem<ode_type, state_type, stepper_type>::
            CalculateMeanField(this->_x.begin()+_node_indices[i],
                               this->_x.begin()+_node_indices[i+1]));
      }
      return mean_field;
    }


    /*!
     * \brief Calculates the coordinates on a sphere of the same dimension as
     * the phase space. If the dimension is 1, the corrdinates will be wrapped
     * around the unit circle. The first coordinate is the radius and the later
     * ones are the phases. Careful: in 3-d this is not the same as spherical
     * coordinates with polar angle and azimuth!
     */
    virtual matrix_type CalculateMeanFieldSpherical() {
      matrix_type mean_field;
      for (size_t i = 0; i < _node_indices.size() - 1; ++i) {
        mean_field.push_back(GenericSystem<ode_type, state_type, stepper_type>::
            CalculateMeanFieldSpherical(this->_x.begin()+_node_indices[i],
                                        this->_x.begin()+_node_indices[i+1]));
      }
      return mean_field;
    }

  protected:
    node_size_type _node_indices;
    node_size_type _node_sizes;

    /*!
     * The network is initialized to a zero state, without initializing the ode.
     *
     * @param node_sizes a vector containing the size of every single node
     * @param dimension the dimensionality of the oscillators
     */
    GenericNetwork(node_size_type node_sizes, unsigned int dimension)
      : GenericSystem<ode_type, state_type, stepper_type>(0, dimension) {
        _node_indices = CalculateNodeIndices(node_sizes);
        _node_sizes = node_sizes;
        this->Resize(CalculateNumberOscillators(_node_sizes));
      }


    node_size_type CalculateNodeIndices(const node_size_type& node_sizes) {
      node_size_type node_indices = {0};
      unsigned int offset = 0;;
      for (size_t i = 0; i < node_sizes.size(); ++i) {
        node_indices.push_back(offset + node_sizes[i]*this->_d);
        offset += node_sizes[i]*this->_d;
      }
      return node_indices;
    }


    unsigned int CalculateNumberOscillators(const node_size_type& node_sizes) {
      unsigned int N = 0;
      for (size_t i = 0; i < node_sizes.size(); ++i) {
        N += node_sizes[i];
      }
      return N;
    }
};
} // nonlinear_systems

#endif
