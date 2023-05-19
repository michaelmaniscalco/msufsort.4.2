#pragma once

#include <include/convertible_to_span_of_bytes.h>
#include <library/system.h>

#include <vector>
#include <span>
#include <cstdint>


namespace maniscalco
{

    //=========================================================================
    template <typename S>
    void make_suffix_array
    (
        std::span<std::byte const>,
        std::span<S>,
        std::uint32_t
    );


    //=========================================================================
    template <typename T, typename S>
    void make_suffix_array
    (
        T source,
        std::span<S> suffixArray,
        std::uint32_t threadCount
    ) requires convertible_to_span_of_bytes<T>
    {
        make_suffix_array(std::as_bytes(std::span(source.begin(), source.size())), suffixArray, threadCount);
    }


    //=========================================================================
    template <typename T, typename S = std::uint32_t>
    std::vector<S> make_suffix_array
    (
        T source,
        std::uint32_t threadCount
    ) requires convertible_to_span_of_bytes<T>
    {
        std::vector<S> suffixArray(source.size() + 1);
        make_suffix_array<T, S>(source, suffixArray, threadCount);
        return suffixArray;
    }

} // namespace maniscalco