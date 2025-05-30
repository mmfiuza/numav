# Main library sources
set(NUMAV_SOURCES numav.cpp)
set(NUMAV_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../include/numav.hpp)

# Create shared library
add_library(numav SHARED ${NUMAV_SOURCES} ${NUMAV_HEADERS})
target_include_directories(numav PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:include>
)