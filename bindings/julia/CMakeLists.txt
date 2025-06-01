# Julia bindings target
add_library(numav_jl SHARED numav_jl.cpp)

# Use imported targets instead of raw variables
target_link_libraries(numav_jl 
    PRIVATE
        numav
        JlCxx::cxxwrap_julia # automatically handles includes
        JlCxx::cxxwrap_julia_stl
        ${Julia_LIBRARY}
)

# Install settings
install(TARGETS numav_jl
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)