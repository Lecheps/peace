cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

include(ExternalProject)

# Adding customized cmake module for building boost
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

project(peace)


    find_package (Python COMPONENTS Interpreter Development REQUIRED) 
    include(BuildBoost) # Custom module

    include_directories(. ${Boost_INCLUDE_DIR} ${Python_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIR})

    # Build and link the peace module
    add_library(peace SHARED peace.cpp RandGenerator.cpp)
    target_link_libraries(peace  ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} -lpthread)
    add_dependencies(peace Boost)
  
    # Tweaks the name of the library to match what Python expects
    set_target_properties(peace PROPERTIES SUFFIX .so)
    set_target_properties(peace PROPERTIES PREFIX "")   
