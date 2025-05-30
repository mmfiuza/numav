include_directories(${CMAKE_SOURCE_DIR}/third-party/cxxwrap/include)

add_library(numav_julia SHARED numav_julia.cpp)
target_link_libraries(numav_julia ${JLCXX_TARGET} ${JLCXX_STL_TARGET} ${Julia_LIBRARY})

install(
  TARGETS 
  numav_julia
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/jlcxx-examples
)
