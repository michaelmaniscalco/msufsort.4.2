set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(msufsort_demo main.cpp)

target_link_libraries(msufsort_demo 
    fmt::fmt
    include 
    msufsort)