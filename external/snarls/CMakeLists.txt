cmake_minimum_required(VERSION 3.16)
project(dinara_snarls)

# Build libhandlegraph and structures as subdirectories.
# Suppress their shared library builds — we only need static.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(libhandlegraph)
add_subdirectory(structures)

# Vendored snarl finder sources (stripped from vg).
add_library(dinara_snarls STATIC
    vendored/three_edge_connected_components.cpp
    vendored/integrated_snarl_finder.cpp
)

target_include_directories(dinara_snarls PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/vendored
    ${CMAKE_CURRENT_SOURCE_DIR}/libhandlegraph/src/include
    ${CMAKE_CURRENT_SOURCE_DIR}/structures/src/include
)

target_link_libraries(dinara_snarls PUBLIC
    handlegraph_static
    structures_static
)

# Match the parent project's C++ standard.
set_target_properties(dinara_snarls PROPERTIES CXX_STANDARD 17)

# Build position-independent code so the static library can be linked
# into shared libraries (e.g. the Python module).
set_target_properties(dinara_snarls PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(handlegraph_static PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
set_target_properties(structures_static PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
