cmake_minimum_required(VERSION 3.8)
project(compdp)

add_subdirectory(./mylib)
add_executable(main main.cpp)
add_executable(baseline baseline.cpp)

target_link_libraries(main simfrontierlib)
target_include_directories(main PRIVATE path/to/boost/)
target_link_libraries(baseline simfrontierlib)