cmake_minimum_required(VERSION 3.9)
project(arbor-rc-exp2syn-spike LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../cmake")
include(UseLibraryPath)

find_package(arbor REQUIRED)
find_package(netcdf REQUIRED)

add_executable(arbor-rc-exp2syn-spike arbor-rc-exp2syn-spike.cpp)
target_link_libraries(arbor-rc-exp2syn-spike arbor::arbor netcdf::netcdf)
target_include_directories(arbor-rc-exp2syn-spike PRIVATE ${CMAKE_SOURCE_DIR}/../include)

install(TARGETS arbor-rc-exp2syn-spike DESTINATION bin)

