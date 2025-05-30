# CxxWrap configuration
set(JLCXX_PROJECT_NAME JlCxx)
set(JLCXX_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(JLCXX_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Headers and sources
file(GLOB JLCXX_HEADERS "${JLCXX_INCLUDE_DIR}/*.hpp")
file(GLOB JLCXX_SOURCES "${JLCXX_SOURCE_DIR}/*.cpp")
file(GLOB JLCXX_STL_SOURCES "${JLCXX_SOURCE_DIR}/stl_*.cpp")

# Version detection
file(STRINGS "${JLCXX_INCLUDE_DIR}/jlcxx_config.hpp" jlcxx_version_defines
    REGEX "#define JLCXX_VERSION_(MAJOR|MINOR|PATCH)")
foreach(ver ${jlcxx_version_defines})
    if(ver MATCHES "#define JLCXX_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(JLCXX_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}")
    endif()
endforeach()
set(${JLCXX_PROJECT_NAME}_VERSION "${JLCXX_VERSION_MAJOR}.${JLCXX_VERSION_MINOR}.${JLCXX_VERSION_PATCH}")

# Create libraries
add_library(cxxwrap_julia SHARED ${JLCXX_SOURCES} ${JLCXX_HEADERS})
add_library(cxxwrap_julia_stl SHARED ${JLCXX_STL_SOURCES})

# Link and include settings
target_include_directories(cxxwrap_julia PUBLIC
    ${JLCXX_INCLUDE_DIR}
    ${Julia_INCLUDE_DIRS}
)
target_link_libraries(cxxwrap_julia ${Julia_LIBRARY})
target_link_libraries(cxxwrap_julia_stl cxxwrap_julia)

# Library properties
set_target_properties(cxxwrap_julia PROPERTIES
    PUBLIC_HEADER "${JLCXX_HEADERS}"
    COMPILE_DEFINITIONS "JLCXX_EXPORTS"
    VERSION ${${JLCXX_PROJECT_NAME}_VERSION}
    SOVERSION ${JLCXX_VERSION_MAJOR}
)
target_compile_definitions(cxxwrap_julia PUBLIC "JULIA_ENABLE_THREADING")
target_compile_features(cxxwrap_julia PUBLIC cxx_std_20)
target_compile_features(cxxwrap_julia_stl PUBLIC cxx_std_20)

# Installation
install(TARGETS cxxwrap_julia cxxwrap_julia_stl
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/jlcxx
)