#pragma once

#include <span>
#include <cstdint>
#include <type_traits>
#include <concepts>


namespace maniscalco
{

    template <typename T> 
    concept convertible_to_span_of_bytes = requires(T a)
            {
                a[0];
                std::as_bytes(std::span(a.begin(), a.size()));
            } && (sizeof(typename T::value_type) == sizeof(std::byte));

} // namespace maniscalco
