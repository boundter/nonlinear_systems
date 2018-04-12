if(NL_BUILD_TESTS)
  
  set(sources 
    test_generic_system.cpp 
    test_rayleigh_mean_field.cpp
    test_observers.cpp
    test_generic_network.cpp)
  
  
  foreach (sourcefile ${sources})
    string(REPLACE ".cpp" "" sourcename ${sourcefile})
    add_executable(${sourcename} ${sourcefile})
    target_link_libraries(${sourcename} nonlinear_systems)
    add_test(${sourcename} ${sourcename})
  endforeach(sourcefile ${sources}) 

endif()
