cmake_minimum_required(VERSION 4.3.0)
project(tblaze_2b LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 5.0 REQUIRED NO_MODULE)
find_package(OpenMP REQUIRED)

add_library(tb src/model.cpp src/select_k.cpp src/dos.cpp src/read_and_write.cpp)
target_include_directories(tb PUBLIC include)
target_link_libraries(tb PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

add_executable(bands_along_path apps/bands_along_path.cpp)
target_link_libraries(bands_along_path PRIVATE tb)

enable_testing()
add_executable(tests tests/test_main.cpp tests/test_graphene.cpp)
target_include_directories(tests PRIVATE third_party)
target_link_libraries(tests PRIVATE tb)
add_test(NAME tests COMMAND tests)
