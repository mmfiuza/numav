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