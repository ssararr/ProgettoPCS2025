# Headers
list(APPEND src_headers
    ${CMAKE_CURRENT_SOURCE_DIR}/PolyhedraMesh.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Duale.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/triangolazione.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/output.hpp
)

# Sources
list(APPEND src_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Duale.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/triangolazione.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/output.cpp
)


list(APPEND src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(src_sources ${src_sources} PARENT_SCOPE)
set(src_headers ${src_headers} PARENT_SCOPE)
set(src_includes ${src_includes} PARENT_SCOPE)
