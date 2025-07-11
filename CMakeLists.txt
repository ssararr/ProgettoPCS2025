cmake_minimum_required(VERSION 3.9)

project(ProgettoPCS2025 LANGUAGES CXX C VERSION 1.0.0)

# Set C++ Standard
################################################################################
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

# Add dependencies
################################################################################
## Eigen3
find_package(Eigen3 CONFIG REQUIRED)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PUBLIC Eigen3::Eigen)

## Threads
find_package(Threads REQUIRED)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest - MODIFICATO: Scommentato e corretto
find_package(GTest REQUIRED)
list(APPEND ProgettoPCS2025_LINKED_LIBRARIES PRIVATE GTest::GTest GTest::Main)

# Insert Sources
################################################################################
add_subdirectory(src)
list(APPEND ProgettoPCS2025_SOURCES ${src_sources})
list(APPEND ProgettoPCS2025_HEADERS ${src_headers})
list(APPEND ProgettoPCS2025_INCLUDE ${src_includes})

add_subdirectory(ExportParaview)

# Aggiungi la directory dei test - MODIFICATO: Scommentato
add_subdirectory(src_test)
list(APPEND ProgettoPCS2025_TEST_SOURCES ${src_test_sources})
list(APPEND ProgettoPCS2025_TEST_HEADERS ${src_test_headers})
list(APPEND ProgettoPCS2025_TEST_INCLUDE ${src_test_includes})

# Create main executable
################################################################################
add_executable(${PROJECT_NAME}
    main.cpp
    ${ProgettoPCS2025_SOURCES}
    ${ProgettoPCS2025_HEADERS})

target_link_libraries(${PROJECT_NAME} 
    PRIVATE 
    ${ProgettoPCS2025_LINKED_LIBRARIES} 
    ExportParaview)
    
target_include_directories(${PROJECT_NAME} 
    PRIVATE 
    ${ProgettoPCS2025_INCLUDE})
    
target_compile_options(${PROJECT_NAME} 
    PUBLIC 
    -fPIC)

# Create test executable - MODIFICATO: Aggiunto linking esplicito a GTest
add_executable(${PROJECT_NAME}_TEST
    main_test.cpp
    ${ProgettoPCS2025_SOURCES}
    ${ProgettoPCS2025_HEADERS}
    ${ProgettoPCS2025_TEST_SOURCES}
    ${ProgettoPCS2025_TEST_HEADERS})

target_link_libraries(${PROJECT_NAME}_TEST 
    PRIVATE 
    ${ProgettoPCS2025_LINKED_LIBRARIES}
    GTest::GTest 
    GTest::Main
    ExportParaview)
    
target_include_directories(${PROJECT_NAME}_TEST 
    PRIVATE 
    ${ProgettoPCS2025_INCLUDE} 
    ${ProgettoPCS2025_TEST_INCLUDE})
    
target_compile_options(${PROJECT_NAME}_TEST 
    PUBLIC 
    -fPIC)
