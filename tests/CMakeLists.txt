# Auto-discover tests in this directory.
# - Each .cpp file becomes a test executable and CTest entry.
# - Each .py file becomes a CTest entry if Python is available.

file(GLOB TEST_CPP_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")
file(GLOB TEST_PY_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.py")

if(TEST_CPP_SOURCES)
    foreach(test_src IN LISTS TEST_CPP_SOURCES)
        get_filename_component(test_name "${test_src}" NAME_WE)

        add_executable("${test_name}" "${test_src}")
        target_link_libraries("${test_name}" PRIVATE bit_info)

        add_test(NAME "${test_name}" COMMAND "${test_name}")
    endforeach()
endif()

if(TEST_PY_SOURCES)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        foreach(test_py IN LISTS TEST_PY_SOURCES)
            get_filename_component(test_name "${test_py}" NAME_WE)
            add_test(NAME "py_${test_name}" COMMAND "${Python3_EXECUTABLE}" "${test_py}")
        endforeach()
    else()
        message(WARNING "Python tests found, but no Python3 interpreter was detected by CMake.")
    endif()
endif()

if(NOT TEST_CPP_SOURCES AND NOT TEST_PY_SOURCES)
    message(STATUS "No tests found in ${CMAKE_CURRENT_SOURCE_DIR} (add .cpp or .py files).")
endif()
