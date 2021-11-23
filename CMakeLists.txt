cmake_minimum_required(VERSION 3.16.3)
set(CMAKE_CXX_STANDARD 20)

set(PROJECT_NAME msusort)
project(${PROJECT_NAME})

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("**** Build type not set.  defaulting to Release ****")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(
        -g
        -O0
        -march=native
    )
else()
    add_compile_options(
        -O3
        -march=native
    )
endif()


option(MSUFSORT_BUILD_DEMO "Build the CLI demo" OFF)

###############################################################################
# include dependency repos 
###############################################################################

include(FetchContent)

# work-contract
FetchContent_Declare(
    work_contract
    GIT_REPOSITORY https://github.com/michaelmaniscalco/work_contract.git
    GIT_TAG master
    SOURCE_DIR        "${CMAKE_BINARY_DIR}/work_contract-src"
    BINARY_DIR        "${CMAKE_BINARY_DIR}/work_contract-build"
    INSTALL_DIR       "${CMAKE_BINARY_DIR}"
    INSTALL_COMMAND   ""
    )
FetchContent_MakeAvailable(work_contract)
FetchContent_GetProperties(work_contract)


# work-contract
FetchContent_Declare(
    common
    GIT_REPOSITORY https://github.com/michaelmaniscalco/common.git
    GIT_TAG master
    SOURCE_DIR        "${CMAKE_BINARY_DIR}/common-src"
    BINARY_DIR        "${CMAKE_BINARY_DIR}/common-build"
    INSTALL_DIR       "${CMAKE_BINARY_DIR}"
    INSTALL_COMMAND   ""
    )
FetchContent_MakeAvailable(common)
FetchContent_GetProperties(common)

# ranges
FetchContent_Declare(
    range_v3
    GIT_REPOSITORY https://github.com/ericniebler/range-v3.git
    GIT_TAG master
    SOURCE_DIR        "${CMAKE_BINARY_DIR}/range_v3-src"
    BINARY_DIR        "${CMAKE_BINARY_DIR}/range_v3-build"
    INSTALL_DIR       "${CMAKE_BINARY_DIR}"
    INSTALL_COMMAND   ""
)
FetchContent_MakeAvailable(range_v3)
FetchContent_GetProperties(range_v3)
set(_range_v3_include_dir "${CMAKE_BINARY_DIR}/range_v3-src/include/")


# fmt
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG master
    SOURCE_DIR        "${CMAKE_BINARY_DIR}/fmt-src"
    BINARY_DIR        "${CMAKE_BINARY_DIR}/fmt-build"
    INSTALL_DIR       "${CMAKE_BINARY_DIR}"
    INSTALL_COMMAND   ""
)
FetchContent_MakeAvailable(fmt)
FetchContent_GetProperties(fmt)
set(_fmt_include_dir "${CMAKE_BINARY_DIR}/fmt-src/include/")


# type_safe
FetchContent_Declare(
    type_safe
    GIT_REPOSITORY https://github.com/foonathan/type_safe.git
    GIT_TAG main
    SOURCE_DIR        "${CMAKE_BINARY_DIR}/type_safe-src"
    BINARY_DIR        "${CMAKE_BINARY_DIR}/type_safe-build"
    INSTALL_DIR       "${CMAKE_BINARY_DIR}"
    INSTALL_COMMAND   ""
)
FetchContent_MakeAvailable(type_safe)
FetchContent_GetProperties(type_safe)
set(_type_safe_include_dir "${CMAKE_BINARY_DIR}/type_safe-src/include/")


add_subdirectory(src)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

