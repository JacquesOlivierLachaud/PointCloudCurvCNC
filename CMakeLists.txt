cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
Project(PointCloudCurvCNC 
        VERSION 1.0
        DESCRIPTION "Lightweight Curvature Estimation on Point Clouds with Randomized Corrected Curvature Measures"
        LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/recipes)

include(config)

add_executable(curvatures curvatures.cpp)
target_link_libraries(curvatures Eigen3::Eigen polyscope)
target_include_directories(curvatures PRIVATE "${CMAKE_SOURCE_DIR}/deps/")

