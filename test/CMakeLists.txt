if(NL_BUILD_TESTS)

  set(CMAKE_BUILD_TYPE Debug)
  
  set(sources
    # misc
    test_coordinates.cpp
    test_statistics.cpp
    test_helper.cpp
    test_derivative.cpp
    # observers
    test_position_observer.cpp
    test_statistics_observer.cpp
    # odes
    test_harmonic_oscillator_ode.cpp
    test_kuramoto_sakaguchi_ode.cpp
    test_rayleigh_mean_field_ode.cpp
    test_m_kuramoto_sakaguchi_ode.cpp
    test_reduced_m_kuramoto_sakaguchi_watanabe_strogatz_ode.cpp
    # frameworks
    test_generic_system.cpp 
    test_generic_network.cpp
    # systems
    test_m_kuramoto_sakaguchi_system.cpp
    test_rayleigh_mean_field_system.cpp
    # dependent stuff
    test_average_frequency_observer.cpp
    test_order_parameter_observer.cpp
    test_instantaneous_frequency_observer.cpp
    test_minimal_order_parameter_observer.cpp)
  
  
  foreach (sourcefile ${sources})
    string(REPLACE ".cpp" "" sourcename ${sourcefile})
    add_executable(${sourcename} ${sourcefile})
    target_link_libraries(${sourcename} nonlinear_systems)
    add_test(${sourcename} ${sourcename})
  endforeach(sourcefile ${sources}) 

  # special compilation
  add_executable(test_parse_arguments test_parse_arguments.cpp)
  target_link_libraries(test_parse_arguments nonlinear_systems ${Boost_LIBRARIES})
  add_test(test_parse_arguments test_parse_arguments)

  add_executable(test_reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system
                 test_reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system.cpp) 
  target_link_libraries(test_reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system
                        nonlinear_systems ${NLOPT_LIBRARIES})
  add_test(test_reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system
           test_reduced_m_kuramoto_sakaguchi_watanabe_strogatz_system)

endif()
