#pragma once

#include <library/system.h>

#include <numeric>
#include <mutex>
#include <span>
#include <condition_variable>


namespace maniscalco
{

    using suffix_value = std::uint64_t;
    using symbol = std::byte;

    class suffix_state;


    template <typename T>
    class msufsort
    {
    public:

        using suffix_index = T;

        struct configuration
        {
            std::uint32_t threadCount_;
        };

        msufsort
        (
            configuration const &
        );

        void suffix_array
        (
            std::span<symbol const>,
            std::span<suffix_index>
        );

    private:

        static std::int32_t constexpr starting_min_match_length_for_tandem_repeats = (2 + sizeof(suffix_value) + sizeof(suffix_value));
        static auto constexpr min_match_length_before_induced_sort = 26;//18;
        static auto constexpr induced_sort_threshold = (2 + (sizeof(suffix_value) * (min_match_length_before_induced_sort / sizeof(suffix_value))));
        static_assert(((induced_sort_threshold - 2) % sizeof(suffix_value)) == 0);
        static_assert(induced_sort_threshold >= (2 + sizeof(suffix_value)));

        static auto constexpr insertion_sort_threshold = 64;

        static auto constexpr bits_per_suffix_index = (8 * sizeof(suffix_index));
        static suffix_index constexpr induced_sort_flag = (1ull << (bits_per_suffix_index - 1));
        static suffix_index constexpr tandem_repeat_flag = induced_sort_flag;  
        static suffix_index constexpr isa_index_mask = ~tandem_repeat_flag;

        struct stack_frame;

        struct tandem_repeat_info
        {
            tandem_repeat_info
            (
                std::span<suffix_index> range,
                std::int32_t numTerminators,
                std::int32_t tandemRepeatLength
            ):
                range_(range),
                numTerminators_(numTerminators),
                tandemRepeatLength_(tandemRepeatLength)
            {
            }

            std::span<suffix_index> range_;
            std::int32_t            numTerminators_;
            std::int32_t            tandemRepeatLength_;
        };
        
        std::vector<tandem_repeat_info> tandemRepeatStack_;

        std::pair<std::uint32_t, std::uint32_t> partition_tandem_repeats
        (
            std::span<suffix_index>,
            std::int32_t,
            std::vector<tandem_repeat_info> &
        );

        bool has_potential_tandem_repeats
        (
            suffix_value,
            std::array<suffix_value, 2>
        ) const;

        void complete_tandem_repeats
        (
            std::span<symbol const>,
            std::vector<tandem_repeat_info> &
        );

        void complete_tandem_repeat
        (
            std::span<symbol const>,
            std::span<suffix_index>,
            std::int32_t,
            std::int32_t
        );

        std::tuple<suffix_index *, suffix_value, suffix_index *, suffix_value, suffix_index *, suffix_value> select_pivots
        (
            symbol const *,
            std::span<suffix_index>
        );

        void start_async_task
        (
            std::size_t,
            std::function<void()>
        );

        suffix_value get_value
        (
            symbol const *,
            suffix_index
        ) const;

         suffix_value __attribute__ ((noinline)) get_value_over
        (
            symbol const *,
            suffix_index
        ) const;

        void count_suffix_types
        (
            std::span<symbol const>,
            std::span<symbol const>,
            std::array<std::uint32_t, 0x40000> &
        );

        void initial_radix_sort
        (
            std::span<symbol const>,
            std::span<suffix_index>,
            std::array<std::uint32_t, 0x40000> &,
            std::uint32_t
        );

        void initial_radix_sort
        (
            std::span<symbol const>,
            std::span<symbol const>,
            std::span<suffix_index>,
            std::array<std::uint32_t, 0x40000> &
        );

        void multikey_quicksort
        (
            std::span<symbol const>,
            std::span<suffix_index>,
            std::uint32_t,
            std::vector<stack_frame> &,
            std::vector<tandem_repeat_info> &
        );
        
        suffix_state update_suffix_state
        (
            std::span<symbol const>,
            std::span<suffix_index>,
            suffix_state,
            suffix_value,
            std::uint32_t
        );

        void multikey_insertion_sort
        (
            std::span<symbol const>,
            std::span<suffix_index>,
            std::uint32_t,
            std::uint32_t,
            suffix_state,
            suffix_value,
            std::array<suffix_value, 2>,
            std::vector<tandem_repeat_info> &
        );

        void quicksort
        (
            std::span<suffix_index>
        );

        void insertion_sort
        (
            suffix_index *,
            suffix_index *
        );

        void complete_isa
        (
            std::span<suffix_index>
        );

        void complete_induced_sort
        (
            std::span<symbol const>,
            std::span<suffix_index>
        );

        void complete_induced_sort
        (
            std::span<suffix_index>
        );

        bool compare_suffixes
        (
            std::span<symbol const>,
            symbol const *,
            suffix_index indexA,
            suffix_index indexB
        ) const;

        void wait_for_all_async_tasks() const;

        bool single_threaded() const;

        void write_isa
        (
            suffix_index,
            suffix_index
        );

        void mark_as_induced_sort
        (
            suffix_index &
        );

        suffix_index clear_induced_sort
        (
            suffix_index &
        );

        bool is_induced_sort
        (
            suffix_index const &
        ) const;

        suffix_index read_isa
        (
            suffix_index
        ) const;

        void mark_suffixes_for_induced_sort
        (
            std::span<suffix_index>,
            std::uint32_t
        );

        std::uint32_t                                   threadCount_;

        system::work_contract_group                     workContractGroup_;

        std::vector<system::work_contract>              workContracts_;

        std::atomic<std::size_t>                        activeWorkContractCount_{0};

        system::thread_pool                             threadPool_;

        std::condition_variable mutable                 conditionVariable_;

        std::mutex mutable                              mutex_;

        symbol const *                                  getValueEnd_;

        suffix_index *                                  sa_;

        suffix_index *                                  isa_;

    }; // class msufsort

} // namespace maniscalco
