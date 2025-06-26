set(EXPORTPARAVIEW_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.cpp)
set(EXPORTPARAVIEW_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/UCDUtilities.hpp)

# Crea la libreria
add_library(ExportParaview STATIC ${EXPORTPARAVIEW_SOURCES} ${EXPORTPARAVIEW_HEADERS})

# Imposta gli include path
target_include_directories(ExportParaview PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ExportParaview PUBLIC Eigen3::Eigen)
target_include_directories(ExportParaview PUBLIC ${EIGEN3_INCLUDE_DIR})

# Esporta le variabili al parent
set(src_paraview_sources ${EXPORTPARAVIEW_SOURCES} PARENT_SCOPE)
set(src_paraview_headers ${EXPORTPARAVIEW_HEADERS} PARENT_SCOPE)
set(src_paraview_includes ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)