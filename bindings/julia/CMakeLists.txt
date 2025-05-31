# Julia bindings target
add_library(numav_julia SHARED numav_julia.cpp)

# Link dependencies
target_link_libraries(numav_julia 
    numav 
    cxxwrap_julia 
    cxxwrap_julia_stl 
    ${Julia_LIBRARY}
)

# Windows-specific link flags
if(WIN32 OR CMAKE_CROSSCOMPILING)
    target_link_options(numav_julia PRIVATE -Wl,--export-all-symbols)
endif()

# Install settings
install(TARGETS numav_julia
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)