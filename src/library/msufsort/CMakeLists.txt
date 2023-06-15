add_library(msufsort
    ./private/msufsort.cpp
    msufsort.cpp
)


target_link_libraries(msufsort 
    system
    include)


target_include_directories(msufsort
    PUBLIC
        "${_range_v3_include_dir}"
        "${_system_include_dir}"
        "${_type_safe_include_dir}"
        $<BUILD_INTERFACE:${_msufsort_include_dir}>
        $<INSTALL_INTERFACE:include/msufsort>
)


target_compile_features(msufsort PUBLIC cxx_std_20)
