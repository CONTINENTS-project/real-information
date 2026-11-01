cmake_minimum_required(VERSION 3.16)
project(real_information_install_check LANGUAGES CXX)

enable_testing()
find_package(real_information 0.1 CONFIG REQUIRED)

# No include paths, C++ standard, or OpenMP options are set by this consumer.
add_executable(check_cpp main.cpp second.cpp)
target_link_libraries(check_cpp PRIVATE real_information::bit_info)
add_test(NAME installed_cpp COMMAND check_cpp "${CMAKE_CURRENT_SOURCE_DIR}/../real_info.yaml")

include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
    enable_language(Fortran)
    get_target_property(real_info_includes real_information::bit_info INTERFACE_INCLUDE_DIRECTORIES)
    find_file(real_info_fortran_source bit_info_mod.f90
        PATHS ${real_info_includes} NO_DEFAULT_PATH)
    if(NOT real_info_fortran_source)
        message(FATAL_ERROR "Installed Fortran interface source was not found")
    endif()
    add_executable(check_fortran "${real_info_fortran_source}" ../test_fortran_abi.f90)
    target_link_libraries(check_fortran PRIVATE real_information::bit_info)
    add_test(NAME installed_fortran
        COMMAND check_fortran "${CMAKE_CURRENT_SOURCE_DIR}/../fixed_point.yaml")
else()
    message(STATUS "Fortran compiler not found; installed Fortran check skipped")
endif()
