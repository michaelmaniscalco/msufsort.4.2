#include "./msufsort.h"

#include "./private/msufsort.h"


//=============================================================================
template <typename S>
void maniscalco::make_suffix_array
(
    std::span<std::byte const> source,
    std::span<S> suffixArray,
    std::uint32_t threadCount
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
        std::uint32_t
    );
    
    template void make_suffix_array<std::uint32_t>
    (
        std::span<std::byte const>,
        std::span<std::uint32_t>,
        std::uint32_t
    );
}

