cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(myTesting VERSION 0.1.0) # update here

# Include dir
include_directories(/usr/local/include)

# Src
AUX_SOURCE_DIRECTORY(src SRC_FILES)

# Headers
set(PROJECT_SOURCE_DIR "src")
set(PROJECT_INCLUDE_DIR "include/myTesting_bits") # update here

# Source files
set(SOURCE_FILES
    ${PROJECT_INCLUDE_DIR}/standard.hpp
    ${PROJECT_INCLUDE_DIR}/globalConstants.hpp
    ${PROJECT_INCLUDE_DIR}/test.hpp
    ${PROJECT_INCLUDE_DIR}/LineChart.hpp
    ${PROJECT_INCLUDE_DIR}/Histogram.hpp
    ${PROJECT_INCLUDE_DIR}/testing.hpp
    ${PROJECT_INCLUDE_DIR}/textFunctions.hpp

    ${PROJECT_INCLUDE_DIR}/mathematics.hpp

    ${PROJECT_INCLUDE_DIR}/RealFunction.hpp
    ${PROJECT_INCLUDE_DIR}/RealFunction2D.hpp

    ${PROJECT_SOURCE_DIR}/test.cpp
    ${PROJECT_SOURCE_DIR}/LineChart.cpp
    ${PROJECT_SOURCE_DIR}/Histogram.cpp
    ${PROJECT_SOURCE_DIR}/testing.cpp
    ${PROJECT_SOURCE_DIR}/textFunctions.cpp

    ${PROJECT_SOURCE_DIR}/mathematics.cpp
)

# Set up such that XCode organizes the files correctly
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES})

# Add library
add_library(myTesting SHARED ${SOURCE_FILES}) # update here

# Include directories
target_include_directories(myTesting PRIVATE include/myTesting_bits) # update here 2 times

# Install
install(TARGETS myTesting DESTINATION lib) # update here

# Install the headers
install(FILES include/myTesting DESTINATION include) # update here

# Create base directory
install(DIRECTORY include/myTesting_bits DESTINATION include) # update here
