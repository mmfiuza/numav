add_executable(sim_sala_config_05_PRATO main.cpp)
target_link_libraries(sim_sala_config_05_PRATO PRIVATE numav)
set_target_properties(sim_sala_config_05_PRATO PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/simulations_bin"
)
