# Julia bindings target
add_library(numav_julia SHARED numav_julia.cpp)

# Link dependencies
target_link_libraries(numav_julia 
    numav 
    cxxwrap_julia 
    cxxwrap_julia_stl 
    ${Julia_LIBRARY}
)

# Install settings
install(TARGETS numav_julia
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)

# Suppress volatile keyword warnings for the julia headers
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
  target_compile_options(numav_julia PRIVATE -Wno-volatile)
endif()