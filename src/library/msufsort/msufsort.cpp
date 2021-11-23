#include "./msufsort.h"

#include "./private/msufsort.h"


//=============================================================================
template <typename S>
void maniscalco::make_suffix_array
(
    std::span<std::byte const> source,
    std::span<S> suffixArray,
    thread_count threadCount
)
{
    msufsort<S>({.threadCount_ = threadCount}).suffix_array(source, suffixArray);
}


//=============================================================================
namespace maniscalco
{
    
    template void make_suffix_array<std::uint64_t>
    (
        std::span<std::byte const>,
        std::span<std::uint64_t>,
        thread_count
    );
    
    template void make_suffix_array<std::uint32_t>
    (
        std::span<std::byte const>,
        std::span<std::uint32_t>,
        thread_count
    );
}

