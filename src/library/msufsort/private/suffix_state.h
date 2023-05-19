#pragma once

#include "./msufsort.h"

#include <include/endian.h>

#include <range/v3/view/enumerate.hpp>


namespace maniscalco
{

    #pragma pack(push, 1)
    class suffix_state
    {
    public:

        suffix_state() = default;

        suffix_state(std::uint32_t state):state_(state){}

        void update
        (
            suffix_value,
            symbol const *
        );

        bool can_be_induce_sorted() const {return compareState_ == compare_state::less;}

        std::uint32_t get_induced_sort_offset() const{return comparedCharacterCount_;}

    private:
    
        using state_type = std::uint32_t;

        static state_type constexpr state_mask = 0x20100;
        static state_type constexpr state_bstar = 0x20000;

        enum class compare_state : std::uint32_t
        {
            equal,
            less,
            greater
        };

        std::uint32_t   comparedCharacterCount_{2};

        state_type      state_{0xffull << 9};

        compare_state   compareState_{compare_state::equal};

        std::uint32_t   currentSuffixLength_{2};

    }; // struct suffix_state
    #pragma pack(pop)
    
} // namespace maniscalco


//=============================================================================
inline void maniscalco::suffix_state::update
(
    suffix_value input,
    symbol const * parentSuffixBegin
)
{
    input = endian_swap<std::endian::native, std::endian::big>(input);
    auto inputCurrent = (std::byte const *)&input;
    auto inputEnd = inputCurrent + sizeof(suffix_value);
    auto parentSuffixCurrent =  parentSuffixBegin + comparedCharacterCount_;
    auto prev = (std::byte)(state_ & 0xff);

    while (inputCurrent < inputEnd)
    {
        if ((compareState_ == compare_state::equal) && (*parentSuffixCurrent != *inputCurrent))
        {
            compareState_ = (*inputCurrent < *parentSuffixCurrent) ? compare_state::less : compare_state::greater;
            if (compareState_ == compare_state::less)
            {
                comparedCharacterCount_ += (inputEnd - inputCurrent); // done
                return;
            }
        }

        auto next = *(inputCurrent++);
        if (next != prev)
        {        
            state_ <<= 9;
            state_ |= static_cast<std::uint8_t>(next);
            if ((((state_ >> 9) - state_) & state_mask) == state_bstar)     
            {
                parentSuffixCurrent =  parentSuffixBegin;
                compareState_ = compare_state::equal;
                auto firstCharInSuffix = (std::byte)((state_ >> 18) & 0xff);
                if (*parentSuffixCurrent != firstCharInSuffix)
                    compareState_ = (firstCharInSuffix < *parentSuffixCurrent) ? compare_state::less : compare_state::greater;

                auto end = (parentSuffixCurrent++) + currentSuffixLength_ - 1;
                while ((compareState_ == compare_state::equal) && (parentSuffixCurrent < end))
                {
                    if (*parentSuffixCurrent != prev)
                        compareState_ = (prev < *parentSuffixCurrent) ? compare_state::less : compare_state::greater;
                    ++parentSuffixCurrent;
                }
                parentSuffixCurrent = end;

                if ((*parentSuffixCurrent != next) && (compareState_ == compare_state::equal))
                    compareState_ = (next < *parentSuffixCurrent) ? compare_state::less : compare_state::greater;
                    
                if (compareState_ == compare_state::less)
                {
                    comparedCharacterCount_ = (currentSuffixLength_ + (inputEnd - inputCurrent)); // done
                    return;
                }
                state_ |= (0xffull << 18);
                comparedCharacterCount_ = (parentSuffixCurrent -  parentSuffixBegin);
            }
            prev = next;
            currentSuffixLength_ = 2;
        }
        ++currentSuffixLength_;
        ++comparedCharacterCount_;
        ++parentSuffixCurrent;
    }
    comparedCharacterCount_ = (parentSuffixCurrent -  parentSuffixBegin);
}
