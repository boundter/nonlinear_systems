if(NL_BUILD_TESTS)

  set(sources test_generic_system.cpp)

  add_executable(nonlinear_systems_tests ${sources})
  target_link_libraries(nonlinear_systems_tests nonlinear_systems)
  add_test(test_generic_system nonlinear_systems_tests)


endif()
