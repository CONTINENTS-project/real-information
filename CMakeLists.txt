cmake_minimum_required(VERSION 3.16)

project(real_information VERSION 0.1.0 LANGUAGES CXX)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)
find_package(OpenMP 4.5 REQUIRED COMPONENTS CXX)

add_library(bit_info SHARED src/bit_info.cpp)
add_library(real_information::bit_info ALIAS bit_info)

target_include_directories(bit_info
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_compile_features(bit_info PUBLIC cxx_std_20)
target_link_libraries(bit_info PUBLIC OpenMP::OpenMP_CXX)

set_target_properties(bit_info PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    OUTPUT_NAME bit_info
    PREFIX ""
)

install(TARGETS bit_info
    EXPORT real_informationTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES src/bit_info.h src/bit_info_mod.f90
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Configure the installed wrapper without embedding the absolute install prefix.
# Source-tree use remains available through REAL_INFO_LIB_PATH.
file(RELATIVE_PATH REAL_INFO_LIBRARY_RELATIVE_PATH
    "${CMAKE_INSTALL_PREFIX}"
    "${CMAKE_INSTALL_FULL_LIBDIR}/bit_info${CMAKE_SHARED_LIBRARY_SUFFIX}"
)
configure_file(src/real_info.py "${CMAKE_CURRENT_BINARY_DIR}/real_info.py" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/real_info.py" DESTINATION .)

set(real_information_cmake_dir "${CMAKE_INSTALL_LIBDIR}/cmake/real_information")
configure_package_config_file(cmake/real_informationConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/real_informationConfig.cmake"
    INSTALL_DESTINATION "${real_information_cmake_dir}"
)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/real_informationConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
)
install(EXPORT real_informationTargets
    NAMESPACE real_information::
    DESTINATION "${real_information_cmake_dir}"
)
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/real_informationConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/real_informationConfigVersion.cmake"
    DESTINATION "${real_information_cmake_dir}"
)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
