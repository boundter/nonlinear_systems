if(NL_BUILD_TESTS)

  set(CMAKE_BUILD_TYPE Debug)
  
  set(sources
    test_coordinates.cpp
    test_statistics.cpp
    test_harmonic_oscillator_ode.cpp
    test_kuramoto_sakaguchi_ode.cpp
    test_rayleigh_mean_field_ode.cpp
    test_m_kuramoto_sakaguchi_ode.cpp
    test_generic_system.cpp 
    test_generic_network.cpp
    test_rayleigh_mean_field.cpp
    test_observers.cpp
    test_m_kuramoto.cpp)
  
  
  foreach (sourcefile ${sources})
    string(REPLACE ".cpp" "" sourcename ${sourcefile})
    add_executable(${sourcename} ${sourcefile})
    target_link_libraries(${sourcename} nonlinear_systems)
    add_test(${sourcename} ${sourcename})
  endforeach(sourcefile ${sources}) 

endif()
