# required CMake minimum version
cmake_minimum_required(VERSION 3.5.0)

# project info
project(arbf-test)

# include directories (-I)
include_directories(/opt/local/include)

# link directories (-L)
link_directories(/opt/local/lib)

# find all header files and save them into DIR_INC
#aux_source_directory(arbf-test/argf-test/include DIR_INC)

# find all source files and save them into DIR_SRC
aux_source_directory(arbf-test/arbf-test/src DIR_SRC)

# specify targets
add_executable(arbf-test ${DIR_SRC})

# set build type
#set(CMAKE_BUILD_TYPE Release)

# extra compilation flags
add_definitions("-std=c++11 -Wno-deprecated-register")

target_link_libraries(arbf-test jsoncpp)
