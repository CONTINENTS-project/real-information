cmake_minimum_required(VERSION 3.16)

project(real_information LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CTest)

add_library(bit_info SHARED src/bit_info.cpp)

target_include_directories(bit_info
	PUBLIC
		${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Mirror manual build flags: -g -O3 -std=c++20 -fPIC -fopenmp
target_compile_options(bit_info PRIVATE -g -O3 -fopenmp)
target_link_options(bit_info PRIVATE -fopenmp)

set_target_properties(bit_info PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME bit_info
    PREFIX ""
)

install(TARGETS bit_info
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(FILES src/real_info.py
    DESTINATION .
)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()