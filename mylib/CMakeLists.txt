cmake_minimum_required(VERSION 3.8)
project(simfrontierlib)

set(SRCS graph.cpp)

add_library(simfrontierlib STATIC ${SRCS})

target_compile_features(simfrontierlib PUBLIC cxx_std_11)
