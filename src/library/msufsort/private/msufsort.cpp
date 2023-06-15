



#include "./msufsort.h"

#include "./suffix_state.h"

#include <include/endian.h>

#include <range/v3/view/enumerate.hpp>

#include <iostream>
#include <chrono>
#include <algorithm>
#include <mutex>
#include <span>

#define VERBOSE
#define VALIDATE

#define USE_PREFETCH

static auto constexpr prefetch_offset = 24;

#ifdef USE_PREFETCH
    #define PREFETCH_READ(addr, offset)  __builtin_prefetch((char *)addr + ((sizeof(addr[0]) * offset)), 1, 0)
    #define PREFETCH_WRITE(addr, offset)  __builtin_prefetch(addr + offset, 1, 1)
#else
    #define PREFETCH_READ(...)   
    #define PREFETCH_WRITE(...)   
#endif


//=============================================================================
#pragma pack(push, 1)
template <typename T>
struct maniscalco::msufsort<T>::stack_frame
{
    std::span<suffix_index> range_;
    std::uint32_t           matchLength_;
    std::uint32_t           minMatchLengthForTandemRepeat_;
    suffix_state            suffixState_;
    suffix_value            startingPattern_;
    std::array<suffix_value, 2> endingPattern_;
};
#pragma pack(pop)


//=============================================================================
template <typename T>
maniscalco::msufsort<T>::msufsort
(
    // construct msufsort instance configured with the specified number of threads.
    configuration const & config
):    
    threadCount_(config.threadCount_),
    workContractGroup_(1 << 16),
    workContracts_(threadCount_),
    threadPool_({})
{
    std::vector<system::thread_pool::thread_configuration> threads(config.threadCount_);
    for (auto & thread : threads)
        thread.function_ = [this](auto const & stopToken){while (!stopToken.stop_requested()) workContractGroup_.service_contracts();};
    threadPool_ = system::thread_pool({threads});

    for (auto & workContract : workContracts_)
        workContract = workContractGroup_.create_contract(nullptr);
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::wait_for_all_async_tasks
(
) const
{
    std::unique_lock uniqueLock(mutex_);
    conditionVariable_.wait(uniqueLock, [&](){return (activeWorkContractCount_ == 0);});    
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::start_async_task
(
    std::size_t contractIndex,
    std::function<void()> task
)
{
    workContracts_[contractIndex].update(
            [this, task = std::move(task)]()
            {
                task(); 
                if (--this->activeWorkContractCount_ == 0)
                {
                    std::unique_lock uniqueLock(this->mutex_);
                    this->conditionVariable_.notify_one();
                }
            });
    ++activeWorkContractCount_;
    workContracts_[contractIndex].invoke();
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::suffix_array
(
    std::span<symbol const> source,
    std::span<suffix_index> suffixArray
)
{
    #ifdef VERBOSE
        auto start = std::chrono::system_clock::now();
    #endif

    getValueEnd_ = (source.data() + source.size() - sizeof(suffix_value));
    sa_ = suffixArray.data();
    isa_ = (suffixArray.data() + ((source.size() + 1) / 2));

    // count suffix types for each unique two symbol combination from source and do initial two byte radix sort on the B* suffixes
    std::array<std::uint32_t, 0x40000> counters{};
    initial_radix_sort(std::span(source.begin(), source.begin() + source.size()), suffixArray, counters, threadCount_);

    #ifdef VERBOSE
        auto finish = std::chrono::system_clock::now();
        auto elapsed = finish - start;
        std::cout << "elasped time for count = " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << "ms\n";
        auto start2 = std::chrono::system_clock::now();
    #endif

    // sort each of the B* suffix partitions as set up by the initial radix sort
    auto currentIndex = 0;
    if (single_threaded())
    {
        // single threaded multikey quicksort
        std::vector<stack_frame> stack;
        stack.reserve((1 << 10) * 16);
        tandemRepeatStack_.reserve(1024);
        for (auto i = 0x20000; i < 0x30000; ++i)
        {
            auto size = counters[i] - currentIndex;
            if (size > 0)
                multikey_quicksort(source, {suffixArray.data() + currentIndex, size}, 2, stack, tandemRepeatStack_);
            currentIndex += size;
        }
    }
    else
    {
        // multithreaded multikey quicksort
        std::vector<std::span<suffix_index>> partitions;
        for (auto i = 0x20000; i < 0x30000; ++i)
        {
            auto size = counters[i];
            if (size > 0)
                partitions.push_back(std::span(suffixArray.data() + currentIndex, size));
            currentIndex += size;
        }

        std::sort(partitions.begin(), partitions.end(), [](auto const & a, auto const & b){return (a.size() < b.size());});
        std::atomic<std::int32_t> partitionCount = partitions.size();
        std::vector<std::vector<tandem_repeat_info>> localTandemRepeatStack;
        localTandemRepeatStack.resize(workContracts_.size());

        for (auto && [index, workContract] :  ranges::views::enumerate(workContracts_))
        {
            auto multikeyQuicksortAllPartitionsTask = [&, index]
                    (
                    ) mutable
                    {
                        std::vector<stack_frame> stack;
                        stack.reserve((1 << 10) * 16);
                        std::int32_t partitionIndex = --partitionCount;
                        while (partitionIndex >= 0)
                        {
                            multikey_quicksort(source, partitions[partitionIndex], 2, stack, localTandemRepeatStack[index]);
                            partitionIndex = --partitionCount;
                        }
                    };
            start_async_task(index, multikeyQuicksortAllPartitionsTask);
        }
        wait_for_all_async_tasks();
        for (auto const & stack : localTandemRepeatStack)
            if (!stack.empty())
                std::copy(stack.begin(), stack.end(), std::back_inserter(tandemRepeatStack_));
    }

    #ifdef VERBOSE
        auto finish2 = std::chrono::system_clock::now();
        elapsed = finish2 - start2;
        std::cout << "elasped time for quicksort = " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << "ms\n";
        auto start3 = std::chrono::system_clock::now();
    #endif

    std::span bStarSuffixes(suffixArray.begin(), currentIndex); 
    // sort remaining B*s which have been deferred to use induced sorting
    complete_induced_sort(source, bStarSuffixes);
    #ifdef VERBOSE
        auto finish3 = std::chrono::system_clock::now();
        elapsed = finish3 - start3;
        std::cout << "elasped time for induced sort = " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << "ms\n";
        elapsed = finish3 - start;
        std::cout << "total elasped time = " << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() << "ms\n";
    #endif

    #ifdef VALIDATE
        auto errorCount = 0;
        bool update = true;
        auto ptr = suffixArray.data();
        for (auto j = 1; j < currentIndex; ++j)
        {
            update = false;
            if (compare_suffixes(source, source.data(), ptr[j - 1], ptr[j]) == true)
            {
                update = true;
                ++errorCount;
            }
            if (((j & 0xfff) == 1) || (j == (currentIndex - 1)))
                update = true;
            if ((update) && (errorCount > 0))
                std::cout << "[ERRORS: " << errorCount << "]  ";
            if (update)
                std::cout << "verifying " << j << "/" << (currentIndex - 1) << '\r' << std::flush;
        }

        std::sort(suffixArray.begin(), suffixArray.begin() + currentIndex);
        for (auto i = 1; i < currentIndex; ++i)
            errorCount += (suffixArray[i-1] >= suffixArray[i]);

        if (errorCount > 0)
            std::cout << "ERRORS: " << errorCount << "                      \n";
        else
            std::cout << "sort validated                          \n";
    #endif
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::complete_isa
(
    // mark the isa for all completed sorted B*s
    std::span<suffix_index> suffixArray
)
{
    auto begin = suffixArray.data();
    auto current = begin;
    auto end = current + suffixArray.size();

    if (single_threaded())
    {
        while (current < end)
        {
            while ((current < end) && (!is_induced_sort(*current)))
            {
                write_isa(*current, current - sa_);
                ++current;
            }
            while ((++current < end) && (!is_induced_sort(*current)))
                ;
            ++current;
        }
        return;
    }
    auto chunkSize = ((suffixArray.size() + threadCount_) / threadCount_);
    for (auto && [index, workContract] :  ranges::views::enumerate(workContracts_))
    {
        start_async_task(index, [begin, end, chunkSize, this]() mutable
                {
                    // TODO: this should not process the pending induced sorted suffixes (like the single threaded version)
                    end = ((begin + chunkSize) > end) ? end : (begin + chunkSize);
                    while (begin < end)
                    {
                        if (!is_induced_sort(*begin))
                            write_isa(*begin, begin - sa_);
                        ++begin;
                    }
                });
        begin += chunkSize;
    }
    wait_for_all_async_tasks();
}


//=============================================================================
template <typename T>
bool maniscalco::msufsort<T>::single_threaded
(
) const
{
    return (threadCount_ == 1);
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::count_suffix_types
(
    std::span<symbol const> source,
    std::span<symbol const> range,
    std::array<std::uint32_t, 0x40000> & counter
)
{
    if (range.empty())
        return;
    auto current = range.data() + range.size() - 1;
    auto sourceEnd = (source.data() + source.size() - 1);
    auto currentNextIsSentinel = (current == sourceEnd);
    std::uint32_t state = ((currentNextIsSentinel) || (current[0] > current[1]));

    if ((!currentNextIsSentinel) && (current[0] <= current[1]))
    {
        while ((++current < sourceEnd) && (current[0] == current[1]))
            ;
        state |= ((current == sourceEnd) || (current[0] > current[1]));  
        current = range.data() + range.size() - 1; 
        state <<= ((current[0] != current[1]) | ((state & 0x01) == 0));
        state |= (current[0] > current[1]);  
    }

    auto prev0 = 0;
    auto prev1 = 0;
    auto prev2 = 0;

    while (true)
    {
        auto k = ((state & 0x03) << 16) | endian_swap<std::endian::native, std::endian::big>(*(std::uint16_t const *)current);
        PREFETCH_WRITE(counter.data(), k);
        ++counter[prev0];
        prev0 = prev1;
        prev1 = prev2;
        prev2 = k;
        if (--current < range.data())
            break;
        state <<= ((current[0] != current[1]) | ((state & 0x01) == 0));
        state |= (current[0] > current[1]);        
    }
    counter[0] -= 3;
    counter[prev0]++;
    counter[prev1]++;
    counter[prev2]++;
}


//=============================================================================
template <typename T>
inline void maniscalco::msufsort<T>::write_isa
(
    suffix_index saIndex,
    suffix_index value
)
{
    auto n = (saIndex & isa_index_mask) / 2;
    if (n == 31010042)
        int y = 9;
    isa_[(saIndex & isa_index_mask) >> 1] = value;
}


//=============================================================================
template <typename T>
inline auto maniscalco::msufsort<T>::read_isa
(
    suffix_index saIndex
) const -> suffix_index
{
    return (isa_[(saIndex & isa_index_mask) >> 1]);
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::initial_radix_sort
(
    std::span<symbol const> source,
    std::span<suffix_index> suffixArray,
    std::array<std::uint32_t, 0x40000> & counters,
    std::uint32_t threadCount
)
{
    // count suffix types for each unique two symbol combination from source
    if (single_threaded())
    {
        // single threaded count
        count_suffix_types(source, source, counters);
        auto counterIndex = (endian_swap<std::endian::native, std::endian::big>(*(std::uint16_t const *)(source.data() + source.size() - 1)));
        ++counters[counterIndex];

        std::array<std::uint32_t, 4> c{0,0,0,0};
        for (auto i = 0; i < 0x10000; ++i)
        {
            c[0] += counters[i];
            c[1] += counters[i + 0x10000];
            c[2] += counters[i + 0x20000];
            c[3] += counters[i + 0x30000];
        }
        // merge subcounters
        std::uint32_t subtotal = 0;
        for (auto i = 0x20000; i < 0x30000; ++i)
        {
            auto temp = counters[i];
            counters[i] = subtotal;
            subtotal += temp;
        }
        // single threaded distribute
        initial_radix_sort(source, source, suffixArray, counters);
    }
    else
    {
        // multi threaded count
        auto chunkSize = (std::uint32_t)((source.size() + threadCount - 1) / threadCount);
        std::vector<std::array<std::uint32_t, 0x40000>> subCounters(threadCount);
        for (auto && [index, subCounter] : ranges::views::enumerate(subCounters))
        {
            auto begin = source.begin() + (index * chunkSize);
            auto size = ((begin + chunkSize) <= source.end()) ? chunkSize : std::distance(begin, source.end());
            start_async_task(index, [&source, &subCounter, range = std::span(begin, size), this](){this->count_suffix_types(source, range, subCounter);});
        }
        wait_for_all_async_tasks();
        auto counterIndex = endian_swap<std::endian::native, std::endian::big>(*(std::uint16_t const *)(source.data() + source.size() - 1));
        ++subCounters.back()[counterIndex];

        for (auto i = 0x20000; i < 0x30000; ++i)
        {
            for (auto & subCounter : subCounters)
                counters[i] += subCounter[i];
        }
        // merge subcounters
        std::uint32_t subtotal = 0;
        for (auto i = 0x20000; i < 0x30000; ++i)
        {        
            for (auto & subCounter : subCounters)
            {
                auto temp = subCounter[i];
                subCounter[i] = subtotal;
                subtotal += temp;
            }
        }
        // multi threaded initial radix sort
        for (auto && [index, subCounter] : ranges::views::enumerate(subCounters))
        {
            auto begin = source.begin() + (index * chunkSize);
            auto size = ((begin + chunkSize) <= source.end()) ? chunkSize : std::distance(begin, source.end());
            start_async_task(index, [&source, &subCounter, &suffixArray, range = std::span(begin, size), this](){this->initial_radix_sort(source, range, suffixArray, subCounter);});
        }
        wait_for_all_async_tasks();
    }
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::initial_radix_sort
(
    std::span<symbol const> source,
    std::span<symbol const> range,
    std::span<suffix_index> suffixArray,
    std::array<std::uint32_t, 0x40000> & counter
)
{
    if (range.empty())
        return;
    auto current = range.data() + range.size() - 1;
    auto sourceEnd = (source.data() + source.size() - 1);
    auto currentNextIsSentinel = (current == sourceEnd);
    // determine which suffix type for last (right most) suffix within the range
    std::uint32_t state = ((currentNextIsSentinel) || (current[0] > current[1]));

    if ((!currentNextIsSentinel) && (current[0] <= current[1]))
    {
        while ((++current < sourceEnd) && (current[0] == current[1]))
            ;
        state |= ((current == sourceEnd) || (current[0] > current[1]));  
        current = range.data() + range.size() - 1; 
        state <<= ((current[0] != current[1]) | ((state & 0x01) == 0));
        state |= (current[0] > current[1]);
    }
    // iterate over all suffixes within the range and determine which suffix type
    // each is.  iteration moves from last to first (right to left) over the range.
    suffix_index suffixIndex = std::distance(source.data(), current);

    auto prev = 0;
    auto prevSuffixIndex = 0;

    while (true)
    {
        if ((state & 0x03) == 0x02)
        {
            auto currentValue = 0x20000ull + endian_swap<std::endian::native, std::endian::big>(*(std::uint16_t const *)current);
            PREFETCH_WRITE(counter.data(), currentValue);
            if (prev != 0)
                suffixArray[counter[prev]++] = prevSuffixIndex;
            prev = currentValue;
            prevSuffixIndex = suffixIndex;
        }
        if (--current < range.data())
            break;
        --suffixIndex;
        state <<= ((current[0] != current[1]) | ((state & 0x01) == 0));
        state |= (current[0] > current[1]);        
    }
    if (prev != 0)
        suffixArray[counter[prev]++] = prevSuffixIndex;
}


//==============================================================================
template <typename T>
bool maniscalco::msufsort<T>::compare_suffixes
(
    std::span<symbol const> source,
    symbol const * inputBegin,
    suffix_index indexA,
    suffix_index indexB
) const
{
    if (indexA > indexB)
        return !compare_suffixes(source, inputBegin, indexB, indexA);
    auto inputCurrentA = inputBegin + indexA;
    auto inputCurrentB = inputBegin + indexB;
    while ((inputCurrentB <= getValueEnd_) && (*(suffix_value const *)inputCurrentB == *(suffix_value const *)inputCurrentA))
    {
        inputCurrentB += sizeof(suffix_value);
        inputCurrentA += sizeof(suffix_value);
    }
    if (inputCurrentB >= (source.data() + source.size()))
        return true;
    return (get_value(inputCurrentA, 0) >= get_value(inputCurrentB, 0));
}


//==============================================================================
template <typename T>
inline void maniscalco::msufsort<T>::mark_as_induced_sort
(
    suffix_index & suffixIndex
)
{
    suffixIndex |= induced_sort_flag;
}


//==============================================================================
template <typename T>
inline auto maniscalco::msufsort<T>::clear_induced_sort
(
    suffix_index & suffixIndex
) -> suffix_index
{
    return (suffixIndex &= ~induced_sort_flag);
}


//==============================================================================
template <typename T>
inline bool maniscalco::msufsort<T>::is_induced_sort
(
    suffix_index const & suffixIndex
) const
{
    return (suffixIndex & induced_sort_flag);
}


//==============================================================================
template <typename T>
inline void maniscalco::msufsort<T>::mark_suffixes_for_induced_sort
(
    // mark all suffixes within this partition as deferred induced sort
    // mark isa for first suffix with the offset from the start of the suffix
    // to the suffix which will be used to induce the final sorted order
    // of all suffixes within this partition.
    std::span<suffix_index> partition,
    std::uint32_t inducedSortOffset
)
{
    write_isa(partition.front(), inducedSortOffset);
    mark_as_induced_sort(partition.front());
    mark_as_induced_sort(partition.back());
}


//==============================================================================
template <typename T>
inline void maniscalco::msufsort<T>::multikey_insertion_sort
(
    std::span<symbol const> source,
    std::span<suffix_index> range,
    std::uint32_t currentMatchLength,
    std::uint32_t minMatchLengthForTandemRepeat,
    suffix_state suffixState,
    suffix_value startingPattern,
    std::array<suffix_value, 2> endingPattern,
    std::vector<tandem_repeat_info> & tandemRepeatStack
)
{
    std::int32_t partitionSize = range.size();
    if (partitionSize == 2)
    {
        if (compare_suffixes(source, source.data() + currentMatchLength, range[0], range[1]))
            std::swap(range[0], range[1]);
        return;
    }    

    auto partitionBegin = range.data();
    if (currentMatchLength >= minMatchLengthForTandemRepeat)
    {
        if (currentMatchLength == starting_min_match_length_for_tandem_repeats)
            startingPattern = get_value(source.data(), *partitionBegin);
        if (has_potential_tandem_repeats(startingPattern, endingPattern))
        {
            auto [tandemRepeatCount, minTandemRepeatLength] = partition_tandem_repeats(range, currentMatchLength, tandemRepeatStack);
            minMatchLengthForTandemRepeat = minTandemRepeatLength;
            if (tandemRepeatCount > 0)
            {
                partitionSize -= tandemRepeatCount;
                partitionBegin += tandemRepeatCount;
                if (partitionSize <= 2)
                {
                    if (partitionSize == 2)
                        if (compare_suffixes(source, source.data() + currentMatchLength, partitionBegin[0], partitionBegin[1]))
                            std::swap(partitionBegin[0], partitionBegin[1]);
                    return;
                }    
            }
        }
    }

    suffix_value cachedValue[insertion_sort_threshold];
    while (true)
    {
        auto sourceOffsetByMatchLength = source.data() + currentMatchLength;
        cachedValue[0] = get_value(sourceOffsetByMatchLength, partitionBegin[0]);
        for (std::int32_t i = 1; i < partitionSize; ++i)
        {
            auto currentIndex = partitionBegin[i];
            suffix_value currentValue = get_value(sourceOffsetByMatchLength, currentIndex);
            auto j = i;
            while ((j > 0) && (cachedValue[j - 1] > currentValue))
            {
                cachedValue[j] = cachedValue[j - 1];
                partitionBegin[j] = partitionBegin[j - 1];
                --j;
            }
            cachedValue[j] = currentValue;
            partitionBegin[j] = currentIndex;
        }

        currentMatchLength += (std::int32_t)sizeof(suffix_value);
        if (cachedValue[0] != cachedValue[partitionSize - 1])
            break;
        // this part prevents stack overflow in the case of all suffixes having very long common match length
        suffixState = update_suffix_state(source, {partitionBegin, (unsigned)partitionSize}, suffixState, cachedValue[0], currentMatchLength);
        if (suffixState.can_be_induce_sorted())
            return;
    }

    for (auto i = 0; i < partitionSize; )
    {
        auto start = i;
        auto startValue = cachedValue[i++];
        while ((i < partitionSize) && (cachedValue[i] == startValue))
            ++i;
        auto n = (i - start);
        if (n > 1)
        {
            auto subPartition = std::span(partitionBegin + start, n);
            auto nextSuffixState = update_suffix_state(source, subPartition, suffixState, startValue, currentMatchLength);
            if (!nextSuffixState.can_be_induce_sorted())
               multikey_insertion_sort(source, subPartition, currentMatchLength, minMatchLengthForTandemRepeat, nextSuffixState, startingPattern, {endingPattern[0], startValue}, tandemRepeatStack);
        }
    }
}


//=============================================================================
template <typename T>
inline auto maniscalco::msufsort<T>::update_suffix_state
(
    std::span<symbol const> source,
    std::span<suffix_index> suffixArray,
    suffix_state suffixState,
    suffix_value suffixValue,
    std::uint32_t suffixLength
) -> suffix_state
{
    if (suffixLength >= induced_sort_threshold)
    {
        if (suffixLength > induced_sort_threshold)
        {
            // ongoing state
            suffixState.update(suffixValue, source.data() + suffixArray[0]);
            if (suffixState.can_be_induce_sorted())
                mark_suffixes_for_induced_sort(suffixArray, suffixLength - suffixState.get_induced_sort_offset());
            return suffixState;
        }
        // initial state
        auto suffixIndex = suffixArray[0];
        suffixState = {(0xff << 9) | (std::uint32_t)source[suffixIndex + 1]};
        auto suffixBegin = source.data() + suffixIndex;
        auto suffixEnd = suffixBegin + suffixLength;
        auto suffixCurrent = suffixBegin + 2;
        while (suffixCurrent < suffixEnd)
        {
            auto suffixValue = get_value(suffixCurrent, 0); 
            suffixCurrent += sizeof(suffix_value);
            suffixState.update(suffixValue, suffixBegin);
            if (suffixState.can_be_induce_sorted())
            {
                auto n = std::distance(suffixBegin, suffixCurrent);
                mark_suffixes_for_induced_sort(suffixArray, n - suffixState.get_induced_sort_offset());
                return suffixState;
            }
        }
    }
    return suffixState;
}


//=============================================================================
template <typename T>
inline bool maniscalco::msufsort<T>::has_potential_tandem_repeats
(
    suffix_value startingPattern,
    std::array<suffix_value, 2> endingPattern
) const
{
    std::int8_t const * end = (std::int8_t const *)endingPattern.data();
    std::int8_t const * begin = end + sizeof(suffix_value);
    while (begin > end)
        if (*(suffix_value const *)--begin == *(suffix_value *)&startingPattern)
            return true;
    return false;
}


//=============================================================================
template <typename T>
std::pair<std::uint32_t, std::uint32_t> maniscalco::msufsort<T>::partition_tandem_repeats
(
    // private:
    // the tandem repeat sort.  determines if the suffixes provided are tandem repeats
    // of other suffixes from within the same group.  If so, sorts the non tandem
    // repeat suffixes and then induces the sorted order of the suffixes which are
    // tandem repeats.
    std::span<suffix_index> range,
    std::int32_t currentMatchLength,
    std::vector<tandem_repeat_info> & tandemRepeatStack
)
{
    auto minDistance = std::numeric_limits<std::int32_t>::max();
    auto partitionBegin = range.data();
    auto partitionEnd = partitionBegin + range.size();
    auto partitionSize = range.size();
    std::sort(partitionBegin, partitionEnd, [](std::int32_t a, std::int32_t b) -> bool{return (a < b);});
    std::int32_t tandemRepeatLength = 0;
    auto const halfCurrentMatchLength = (currentMatchLength >> 1);
    // determine if there are tandem repeats and, if so, what the tandem repeat length is.
    auto previousSuffixIndex = partitionBegin[0];
    for (auto cur = partitionBegin + 1; ((tandemRepeatLength == 0) && (cur < partitionEnd)); ++cur)
    {
        auto currentSuffixIndex = *cur;
        auto distance = currentSuffixIndex - previousSuffixIndex;
        if (distance < minDistance)
            minDistance = distance;
        if ((tandemRepeatLength == 0) && ((previousSuffixIndex + halfCurrentMatchLength) >= currentSuffixIndex))
            tandemRepeatLength = (currentSuffixIndex - previousSuffixIndex);
        previousSuffixIndex = currentSuffixIndex;
    }

    if (tandemRepeatLength == 0)
        return {0, minDistance * 2}; // no tandem repeats were found

    // tandem repeats detected.
    auto terminatorsEnd = partitionEnd - 1;
    previousSuffixIndex = partitionEnd[-1];
    for (auto cur = partitionEnd - 2; cur >= partitionBegin; --cur)
    {
	    auto currentSuffixIndex = *cur;
	    if ((previousSuffixIndex - currentSuffixIndex) == tandemRepeatLength)
            std::swap(*terminatorsEnd--, *cur);// suffix is a tandem repeat
	    previousSuffixIndex = currentSuffixIndex;
    }
    auto numTerminators = (std::distance(partitionBegin, terminatorsEnd) + 1);
    std::reverse(partitionBegin, partitionEnd);
    tandemRepeatStack.push_back(tandem_repeat_info({partitionBegin, (unsigned)std::distance(partitionBegin, partitionEnd)}, 
            (std::int32_t)numTerminators, tandemRepeatLength));
    return {(partitionSize - numTerminators), minDistance * 2};
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::complete_tandem_repeats
(
    std::span<symbol const> source,
    std::vector<tandem_repeat_info> & tandemRepeatStack
)
{
    while (!tandemRepeatStack.empty())
    {
        tandem_repeat_info tandemRepeat = tandemRepeatStack.back();
        tandemRepeatStack.pop_back();
        complete_tandem_repeat(source, tandemRepeat.range_, tandemRepeat.numTerminators_, tandemRepeat.tandemRepeatLength_);
    }
}


//=============================================================================
template <typename T>
inline void maniscalco::msufsort<T>::complete_tandem_repeat
(
    std::span<symbol const> source,
    std::span<suffix_index> range,
    std::int32_t numTerminators,
    std::int32_t tandemRepeatLength
)
{
    auto tandemRepeatFlag = (tandem_repeat_flag | tandemRepeatLength);
    auto partitionBegin = range.data();
    auto partitionEnd = partitionBegin + range.size();
    auto * terminatorsBegin = partitionEnd - numTerminators;

    complete_induced_sort({terminatorsBegin, terminatorsBegin + numTerminators});
    // mark isa for each non terminator as length of tandem repeat
    for (auto cur = partitionBegin; cur < terminatorsBegin; ++cur)
        write_isa(*cur, tandemRepeatFlag);

    // now use sorted order of terminators to determine sorted order of repeats.
    // figure out how many terminators sort before the repeat and how
    // many sort after the repeat.  put them on left and right extremes of the array.
    std::int32_t m = 0;
    std::int32_t a = 0;
    std::int32_t b = numTerminators - 1;
    std::int32_t numTypeA = 0;
    while (a <= b)
    {
        m = (a + b) >> 1;
        if (!compare_suffixes(source, source.data() + terminatorsBegin[m], 0, tandemRepeatLength))
        {
	        numTypeA = m;
	        b = m - 1;
        }
        else
        {
	        numTypeA = m + 1;
	        a = m + 1;
        }
    }
    if (numTypeA > numTerminators)
        numTypeA = numTerminators;
    std::int32_t numTypeB = (numTerminators - numTypeA);
    // move A terminators to front of partition
    for (std::int32_t i = 0; i < numTypeA; ++i)
        partitionBegin[i] = terminatorsBegin[i];

    // type A repeats
    auto current = partitionBegin;
    auto currentEnd = current + numTypeA;
    auto next = currentEnd;
    while (current != currentEnd)
    { 
        while (current != currentEnd)
        {
            auto index = *current++;
            if (index >= tandemRepeatLength)
            {
                auto potentialTandemRepeatIndex = index - tandemRepeatLength;
                auto isaValue = read_isa(potentialTandemRepeatIndex);
                if (isaValue == tandemRepeatFlag)
                    *(next++) = potentialTandemRepeatIndex;         
            }
        }
        currentEnd = next;
    }

    // type B repeats
    current = partitionEnd - 1;
    currentEnd = current - numTypeB;
    next = currentEnd;
    while (current != currentEnd)
    { 
        while (current != currentEnd)
        {
            auto index = *current--;
            if (index >= tandemRepeatLength)
            {
                auto potentialTandemRepeatIndex = index - tandemRepeatLength;
                auto isaValue = read_isa(potentialTandemRepeatIndex);
                if (isaValue == tandemRepeatFlag)
                    *(next--) = potentialTandemRepeatIndex;         
            }
        }
        currentEnd = next;
    }
}

/*
//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::multikey_radix_sort
(
    std::span<symbol const> source,
    std::span<suffix_index> range[2],
    std::uint32_t currentMatchLength
)
{
}
*/


//=============================================================================
template <typename T>
auto maniscalco::msufsort<T>::select_pivots
(
    symbol const * source,
    std::span<suffix_index> partition
) -> std::tuple<suffix_index *, suffix_value, suffix_index *, suffix_value, suffix_index *, suffix_value>
{
    auto stride = (partition.size() >> 2);
    auto pivotCandidate1 = partition.data();
    auto pivotCandidateValue1 = get_value(source, *pivotCandidate1);
    auto pivotCandidate2 = pivotCandidate1 + stride;
    auto pivotCandidateValue2 = get_value(source, *pivotCandidate2);
    auto pivotCandidate3 = pivotCandidate2 + stride;
    auto pivotCandidateValue3 = get_value(source, *pivotCandidate3);
    auto pivotCandidate4 = pivotCandidate3 + stride;
    auto pivotCandidateValue4 = get_value(source, *pivotCandidate4);
    auto pivotCandidate5 = pivotCandidate4 + stride - 1;
    auto pivotCandidateValue5 = get_value(source, *pivotCandidate5);
    if (pivotCandidateValue1 > pivotCandidateValue2)
        std::swap(*pivotCandidate1, *pivotCandidate2), std::swap(pivotCandidateValue1, pivotCandidateValue2);
    if (pivotCandidateValue1 > pivotCandidateValue3)
        std::swap(*pivotCandidate1, *pivotCandidate3), std::swap(pivotCandidateValue1, pivotCandidateValue3);
    if (pivotCandidateValue2 > pivotCandidateValue3)
        std::swap(*pivotCandidate2, *pivotCandidate3), std::swap(pivotCandidateValue2, pivotCandidateValue3);
    if (pivotCandidateValue4 > pivotCandidateValue5)
        std::swap(*pivotCandidate4, *pivotCandidate5), std::swap(pivotCandidateValue4, pivotCandidateValue5);
    if (pivotCandidateValue1 > pivotCandidateValue4)
        std::swap(*pivotCandidate1, *pivotCandidate4), std::swap(pivotCandidateValue1, pivotCandidateValue4);
    if (pivotCandidateValue3 > pivotCandidateValue4)
        std::swap(*pivotCandidate3, *pivotCandidate4), std::swap(pivotCandidateValue3, pivotCandidateValue4);
    if (pivotCandidateValue2 > pivotCandidateValue5)
        std::swap(*pivotCandidate2, *pivotCandidate5), std::swap(pivotCandidateValue2, pivotCandidateValue5);
    if (pivotCandidateValue2 > pivotCandidateValue3)
        std::swap(*pivotCandidate2, *pivotCandidate3), std::swap(pivotCandidateValue2, pivotCandidateValue3);
    if (pivotCandidateValue4 > pivotCandidateValue5)
        std::swap(*pivotCandidate4, *pivotCandidate5), std::swap(pivotCandidateValue4, pivotCandidateValue5);
    return {pivotCandidate1, pivotCandidateValue1, pivotCandidate3, pivotCandidateValue3, pivotCandidate5, pivotCandidateValue5};
}

/*
//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::multikey_quicksort_r
(
    std::span<symbol const> source,
    std::span<suffix_index> range[2],
    std::uint32_t currentMatchLength,
    std::vector<stack_frame> & stack,
    std::vector<tandem_repeat_info> & tandemRepeatStack
)
{
    stack.push_back({range, currentMatchLength, starting_min_match_length_for_tandem_repeats, 0});
    suffix_value startingPattern;
    std::array<suffix_value, 2> endingPattern;

    while (!stack.empty())
    {
        auto & next = stack.back();
        auto currentRange = next.range_;
        auto partitionBegin = currentRange.data();
        auto partitionEnd = partitionBegin + currentRange.size();
        auto suffixState = next.suffixState_;
        auto minMatchLengthForTandemRepeat = next.minMatchLengthForTandemRepeat_;
        currentMatchLength = next.matchLength_;
        endingPattern = next.endingPattern_;
        stack.pop_back();

        std::uint64_t partitionSize = std::distance(partitionBegin, partitionEnd);
        if (currentMatchLength >= minMatchLengthForTandemRepeat)
        {
            if (currentMatchLength == starting_min_match_length_for_tandem_repeats)
                startingPattern = get_value(source.data(), *partitionBegin);
            if ((partitionSize > 1) && (has_potential_tandem_repeats(startingPattern, endingPattern)))
            {
                auto [tandemRepeatCount, minTandemRepeatDistance] = partition_tandem_repeats(currentRange, currentMatchLength, tandemRepeatStack);
                partitionBegin += tandemRepeatCount;
                minMatchLengthForTandemRepeat = minTandemRepeatDistance;
                partitionSize = std::distance(partitionBegin, partitionEnd);
            }
        }
        if (partitionSize < insertion_sort_threshold)
        {
            if (partitionSize > 1)
               multikey_insertion_sort(source, {partitionBegin, partitionSize}, currentMatchLength, minMatchLengthForTandemRepeat,
                        suffixState, startingPattern, endingPattern, tandemRepeatStack);
            continue;
        }

        // select three pivots and partition seven ways
        auto offsetInputBegin = source.data() + currentMatchLength;
        auto [pivotCandidate1, pivot1, pivotCandidate2, pivot2, pivotCandidate3, pivot3] = select_pivots(offsetInputBegin, {partitionBegin, partitionSize});
        // destination partition count
        auto p[4] = [pivot1, pivot2, pivot3, 0xffffffffffffffff];
        std::array<std::uint64_t, 8> bucketCount;
        for (auto i = 0; i < partitionSize; ++i)
        {
            auto c = cached[i];
            auto n = (c > pivot2);
            n <<= 2;
            n |= ((c > p[1 + n]) << 1);
            n |= (c > p[n]);
            ++bucketCount[n];
        }
        auto t = 0;
        for (auto & _ : bucketCount)
        {
            auto temp = _;
            _ = t;
            t += temp;
        }
        // distribute
        for (auto i = 0; i < partitionSize; ++i)
        {
            auto c = cached[i];
            auto n = (c > pivot2);
            n <<= 2;
            n |= ((c > p[1 + n]) << 1);
            n |= (c > p[n]);
            if (n & 0x01)
            {
                // is a pivot match
            }
            else
            {
                
            }
        }        

    }
}
*/

//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::multikey_quicksort
(
    std::span<symbol const> source,
    std::span<suffix_index> range,
    std::uint32_t currentMatchLength,
    std::vector<stack_frame> & stack,
    std::vector<tandem_repeat_info> & tandemRepeatStack
)
{
    stack.push_back({range, currentMatchLength, starting_min_match_length_for_tandem_repeats, 0});
    suffix_value startingPattern;
    std::array<suffix_value, 2> endingPattern;

    while (!stack.empty())
    {
        auto & next = stack.back();
        auto currentRange = next.range_;
        auto partitionBegin = currentRange.data();
        auto partitionEnd = partitionBegin + currentRange.size();
        auto suffixState = next.suffixState_;
        auto minMatchLengthForTandemRepeat = next.minMatchLengthForTandemRepeat_;
        currentMatchLength = next.matchLength_;
        endingPattern = next.endingPattern_;
        stack.pop_back();

        std::uint64_t partitionSize = std::distance(partitionBegin, partitionEnd);
        if (currentMatchLength >= minMatchLengthForTandemRepeat)
        {
            if (currentMatchLength == starting_min_match_length_for_tandem_repeats)
                startingPattern = get_value(source.data(), *partitionBegin);
            if ((partitionSize > 1) && (has_potential_tandem_repeats(startingPattern, endingPattern)))
            {
                auto [tandemRepeatCount, minTandemRepeatDistance] = partition_tandem_repeats(currentRange, currentMatchLength, tandemRepeatStack);
                partitionBegin += tandemRepeatCount;
                minMatchLengthForTandemRepeat = minTandemRepeatDistance;
                partitionSize = std::distance(partitionBegin, partitionEnd);
            }
        }
        if (partitionSize < insertion_sort_threshold)
        {
            if (partitionSize > 1)
               multikey_insertion_sort(source, {partitionBegin, partitionSize}, currentMatchLength, minMatchLengthForTandemRepeat,
                        suffixState, startingPattern, endingPattern, tandemRepeatStack);
            continue;
        }

        // select three pivots and partition seven ways
        auto offsetInputBegin = source.data() + currentMatchLength;
        auto [pivotCandidate1, pivot1, pivotCandidate2, pivot2, pivotCandidate3, pivot3] = select_pivots(offsetInputBegin, {partitionBegin, partitionSize});
        auto curSuffix = partitionBegin;
        auto beginPivot1 = partitionBegin;
        auto endPivot1 = partitionBegin;
        auto beginPivot2 = partitionBegin;
        auto endPivot2 = partitionEnd - 1;
        auto beginPivot3 = endPivot2;
        auto endPivot3 = endPivot2;

        std::swap(*curSuffix++, *pivotCandidate1);
        beginPivot2 += (pivot1 != pivot2);
        endPivot1 += (pivot1 != pivot2);
        std::swap(*curSuffix++, *pivotCandidate2);
        if (pivot2 != pivot3)
        {
            std::swap(*endPivot2--, *pivotCandidate3);
            --beginPivot3;
        }

        auto currentValue = get_value(offsetInputBegin, *curSuffix);
        auto nextValue = get_value(offsetInputBegin, curSuffix[1]);
        auto nextDValue = get_value(offsetInputBegin, *endPivot2);

        while (curSuffix <= endPivot2)
        {
            if (currentValue <= pivot2)
            {
                PREFETCH_READ(offsetInputBegin, curSuffix[prefetch_offset]);
                auto temp = nextValue;
                nextValue = get_value(offsetInputBegin, curSuffix[2]);
                if (currentValue < pivot2)
                {
                    std::swap(*beginPivot2, *curSuffix);
                    if (currentValue <= pivot1)
                    {
                        if (currentValue < pivot1)
                            std::swap(*beginPivot1++, *beginPivot2);
                        std::swap(*endPivot1++, *beginPivot2); 
                    }
                    ++beginPivot2; 
                }
                ++curSuffix;
                currentValue = temp;
            }
            else
            {
                PREFETCH_READ(offsetInputBegin, endPivot2[-prefetch_offset]);
                auto temp = get_value(offsetInputBegin, endPivot2[-1]);
                std::swap(*endPivot2, *curSuffix);
                if (currentValue >= pivot3)
                {
                    if (currentValue > pivot3)
                        std::swap(*endPivot2, *endPivot3--);
                    std::swap(*endPivot2, *beginPivot3--);
                }
                --endPivot2;
                currentValue = nextDValue;
                nextDValue = temp;
            }
        }
        ++endPivot3;
        ++beginPivot3;
        ++endPivot2;
        auto nextMatchLength = (currentMatchLength + (std::uint32_t)sizeof(suffix_value));

        if (auto subPartition = std::span(endPivot3, partitionEnd - endPivot3); subPartition.size() > 1)
            stack.push_back({subPartition, currentMatchLength, minMatchLengthForTandemRepeat, suffixState, startingPattern, endingPattern});

        if (auto subPartition = std::span(beginPivot3, endPivot3 - beginPivot3); subPartition.size() > 1)
            if (auto nextSuffixState = update_suffix_state(source, subPartition, suffixState, pivot3, nextMatchLength); !nextSuffixState.can_be_induce_sorted())
                stack.push_back({subPartition, nextMatchLength, minMatchLengthForTandemRepeat, nextSuffixState, startingPattern, {endingPattern[1], pivot3}});

        if (auto subPartition = std::span(endPivot2, beginPivot3 - endPivot2); subPartition.size() > 1)
            stack.push_back({subPartition, currentMatchLength, minMatchLengthForTandemRepeat, suffixState, startingPattern, endingPattern});

        if (auto subPartition = std::span(beginPivot2, endPivot2 - beginPivot2); subPartition.size() > 1)
            if (auto nextSuffixState = update_suffix_state(source, subPartition, suffixState, pivot2, nextMatchLength); !nextSuffixState.can_be_induce_sorted())
                stack.push_back({subPartition, nextMatchLength, minMatchLengthForTandemRepeat, nextSuffixState, startingPattern, {endingPattern[1], pivot2}});

        if (auto subPartition = std::span(endPivot1, beginPivot2 - endPivot1); subPartition.size() > 1)
            stack.push_back({std::span(endPivot1, beginPivot2 - endPivot1), currentMatchLength, minMatchLengthForTandemRepeat, suffixState, startingPattern, endingPattern});

        if (auto subPartition = std::span(beginPivot1, endPivot1 - beginPivot1); subPartition.size() > 1)
            if (auto nextSuffixState = update_suffix_state(source, subPartition, suffixState, pivot1, nextMatchLength); !nextSuffixState.can_be_induce_sorted())
                stack.push_back({subPartition, nextMatchLength, minMatchLengthForTandemRepeat, nextSuffixState, startingPattern, {endingPattern[1], pivot1}});

        if (auto subPartition = std::span(partitionBegin, beginPivot1 - partitionBegin); subPartition.size() > 1)
            stack.push_back({subPartition, currentMatchLength, minMatchLengthForTandemRepeat, suffixState, startingPattern, endingPattern});
    }
}


//==============================================================================
template <typename T>
void maniscalco::msufsort<T>::complete_induced_sort
(
    std::span<symbol const> source,
    std::span<suffix_index> suffixArray
)
{
    std::sort(tandemRepeatStack_.begin(), tandemRepeatStack_.end(), 
            [](auto const & a, auto const & b)
            {
                if (a.range_.end() != b.range_.end())
                    return (a.range_.end() < b.range_.end());
                return (a.range_.begin() < b.range_.begin());
            });

    auto curSuffixArray = suffixArray.begin();
    for (auto curTandemRepeat = tandemRepeatStack_.begin(); curTandemRepeat < tandemRepeatStack_.end(); )
    {
        if (curTandemRepeat->range_.begin() > curSuffixArray)
        {
            // complete all suffixes sorted before the tandem repeats
            complete_isa({curSuffixArray, curTandemRepeat->range_.begin()});
            complete_induced_sort({curSuffixArray, curTandemRepeat->range_.begin()});
        }
        // complete the tandem repeats
        curSuffixArray = curTandemRepeat->range_.end();
        auto start = curTandemRepeat++;
        while ((curTandemRepeat < tandemRepeatStack_.end()) && (curTandemRepeat->range_.end() == start->range_.end()))
            ++curTandemRepeat;
        for (auto n = curTandemRepeat - 1; n >= start; --n)
            complete_tandem_repeat(source, n->range_, n->numTerminators_, n->tandemRepeatLength_);
        complete_isa(start->range_);
    }
    complete_isa({curSuffixArray, suffixArray.end()});
    complete_induced_sort({curSuffixArray, suffixArray.end()});
}


//==============================================================================
template <typename T>
void maniscalco::msufsort<T>::complete_induced_sort
(
    std::span<suffix_index> suffixArray
)
{
    auto current = suffixArray.data();
    auto end = current + suffixArray.size();

    while (current < end)
    {
        while ((current < end) && (!is_induced_sort(*current)))
            ++current;
        if (current < end)
        {
            auto start = current;
            auto offset = read_isa(*start);
            *start = clear_induced_sort(*start) + offset;
            ++current;
            while ((current < end) &&  (!is_induced_sort(*current)))
            {
                PREFETCH_READ(isa_, current[0] >> 1);
                *current += offset;
                ++current;
            }
            *current = clear_induced_sort(*current) + offset;
            ++current;

            if ((current - start) < insertion_sort_threshold)
                insertion_sort(start, current);
            else
                quicksort({start, current});
            while (start < current)
            {
                *start -= offset;
                write_isa(*start, start - sa_);
                start++;
            }
        }
    }
}


//==============================================================================
template <typename T>
inline void maniscalco::msufsort<T>::insertion_sort
(
    suffix_index * begin,
    suffix_index * end
)
{
    std::int32_t partitionSize = std::distance(begin, end);
    if (partitionSize == 2)
    {
        if (read_isa(begin[0]) > read_isa(begin[1]))
            std::swap(begin[0], begin[1]);
        return;
    }

    suffix_index value[insertion_sort_threshold];
    auto n = 0;
    for (auto cur = begin; cur < end; ++cur)
        value[n++] = read_isa(cur[0]);

    for (std::int32_t i = 1; i < partitionSize; ++i)
    {
        auto currentIndex = begin[i];
        auto currentValue = value[i];
        auto j = i;
        while ((j > 0) && (value[j - 1] > currentValue))
        {
            value[j] = value[j - 1];
            begin[j] = begin[j - 1];
            --j;
        }
        value[j] = currentValue;
        begin[j] = currentIndex;
    }
}


//=============================================================================
template <typename T>
void maniscalco::msufsort<T>::quicksort
(
    std::span<suffix_index> partition
)
{
    struct
    {
        suffix_index * begin_;
        suffix_index * end_;
    } stack[1024];

    auto stackTop = stack;
    *(stackTop++) = {partition.data(), partition.data() + partition.size()};

    while (stackTop != stack)
    {
        auto partition = *--stackTop;
        auto partitionBegin = partition.begin_;
        auto partitionEnd = partition.end_;
        auto partitionSize = std::distance(partitionBegin, partitionEnd);

        if (partitionSize < insertion_sort_threshold)
        {
            insertion_sort(partitionBegin, partitionEnd);
            continue;
        }

        // select three pivots
        auto oneSixthOfPartitionSize = (partitionSize / 6);
        auto pivotCandidate1 = partitionBegin + oneSixthOfPartitionSize;
        auto pivotCandidate2 = pivotCandidate1 + oneSixthOfPartitionSize;
        auto pivotCandidate3 = pivotCandidate2 + oneSixthOfPartitionSize;
        auto pivotCandidate4 = pivotCandidate3 + oneSixthOfPartitionSize;
        auto pivotCandidate5 = pivotCandidate4 + oneSixthOfPartitionSize;
        auto pivotCandidateValue1 = read_isa(*pivotCandidate1);
        auto pivotCandidateValue2 = read_isa(*pivotCandidate2);
        auto pivotCandidateValue3 = read_isa(*pivotCandidate3);
        auto pivotCandidateValue4 = read_isa(*pivotCandidate4);
        auto pivotCandidateValue5 = read_isa(*pivotCandidate5);
        if (pivotCandidateValue1 > pivotCandidateValue2)
            std::swap(*pivotCandidate1, *pivotCandidate2), std::swap(pivotCandidateValue1, pivotCandidateValue2);
        if (pivotCandidateValue4 > pivotCandidateValue5)
            std::swap(*pivotCandidate4, *pivotCandidate5), std::swap(pivotCandidateValue4, pivotCandidateValue5);
        if (pivotCandidateValue1 > pivotCandidateValue3)
            std::swap(*pivotCandidate1, *pivotCandidate3), std::swap(pivotCandidateValue1, pivotCandidateValue3);
        if (pivotCandidateValue2 > pivotCandidateValue3)
            std::swap(*pivotCandidate2, *pivotCandidate3), std::swap(pivotCandidateValue2, pivotCandidateValue3);
        if (pivotCandidateValue1 > pivotCandidateValue4)
            std::swap(*pivotCandidate1, *pivotCandidate4), std::swap(pivotCandidateValue1, pivotCandidateValue4);
        if (pivotCandidateValue3 > pivotCandidateValue4)
            std::swap(*pivotCandidate3, *pivotCandidate4), std::swap(pivotCandidateValue3, pivotCandidateValue4);
        if (pivotCandidateValue2 > pivotCandidateValue5)
            std::swap(*pivotCandidate2, *pivotCandidate5), std::swap(pivotCandidateValue2, pivotCandidateValue5);
        if (pivotCandidateValue2 > pivotCandidateValue3)
            std::swap(*pivotCandidate2, *pivotCandidate3), std::swap(pivotCandidateValue2, pivotCandidateValue3);
        if (pivotCandidateValue4 > pivotCandidateValue5)
            std::swap(*pivotCandidate4, *pivotCandidate5), std::swap(pivotCandidateValue4, pivotCandidateValue5);
        auto pivot1 = pivotCandidateValue1;
        auto pivot2 = pivotCandidateValue3;
        auto pivot3 = pivotCandidateValue5;
        // partition seven ways
        auto curSuffix = partitionBegin;
        auto beginPivot1 = partitionBegin;
        auto endPivot1 = partitionBegin;
        auto beginPivot2 = partitionBegin;
        auto endPivot2 = partitionEnd - 1;
        auto beginPivot3 = endPivot2;
        auto endPivot3 = endPivot2;

        std::swap(*(curSuffix++), *pivotCandidate1);
        beginPivot2 += (pivot1 != pivot2);
        endPivot1 += (pivot1 != pivot2);
        std::swap(*(curSuffix++), *pivotCandidate3);
        if (pivot2 != pivot3)
        {
            std::swap(*(endPivot2--), *pivotCandidate5);
            --beginPivot3;
        }

        auto currentValue = read_isa(curSuffix[0]);
        auto nextValue = read_isa(curSuffix[1]);
        auto nextDValue = read_isa(endPivot2[0]);
        while (curSuffix <= endPivot2)
        {
            if (currentValue <= pivot2)
            {
                auto temp = nextValue;
                nextValue = read_isa(curSuffix[2]);
                if (currentValue < pivot2)
                {
                    std::swap(*beginPivot2, *curSuffix);
                    if (currentValue <= pivot1)
                    {
                        if (currentValue < pivot1)
                            std::swap(*beginPivot1++, *beginPivot2);
                        std::swap(*endPivot1++, *beginPivot2); 
                    }
                    ++beginPivot2; 
                }
                ++curSuffix;
                currentValue = temp;
            }
            else
            {
                auto temp = read_isa(endPivot2[-1]);
                std::swap(*endPivot2, *curSuffix);
                if (currentValue >= pivot3)
                {
                    if (currentValue > pivot3)
                        std::swap(*endPivot2, *endPivot3--);
                    std::swap(*endPivot2, *beginPivot3--);
                }
                --endPivot2;
                currentValue = nextDValue;
                nextDValue = temp;
            }
        }

        ++endPivot3;
        ++beginPivot3;
        ++endPivot2;

        if ((partitionEnd - endPivot3) > 1)
            *stackTop++ = {endPivot3, partitionEnd};

        if ((beginPivot3 - endPivot2) > 1)
            *stackTop++ = {endPivot2, beginPivot3};

        if ((beginPivot2 - endPivot1) > 1)
            *stackTop++ = {endPivot1, beginPivot2};

        if ((beginPivot1 - partitionBegin) > 1)
            *stackTop++ = {partitionBegin, beginPivot1};
    }
}


//==============================================================================
template <typename T>
inline auto maniscalco::msufsort<T>::get_value
(
    symbol const * source,
    suffix_index index
) const -> suffix_value
{
    source += index;
    if (source <= getValueEnd_) [[likely]]
        return endian_swap<std::endian::native, std::endian::big>(*(suffix_value const *)(source));
    else [[unlikely]]
        return get_value_over(source, index);
}


//=============================================================================
template <typename T>
auto __attribute__ ((noinline)) maniscalco::msufsort<T>::get_value_over
(
    symbol const * source,
    suffix_index index
) const -> suffix_value
{
    auto over = std::distance(getValueEnd_, source);
    return (over >= sizeof(suffix_value)) ? 0 : (endian_swap<std::endian::native, std::endian::big>(*(suffix_value const *)(getValueEnd_)) << (over * 8));
}


//=============================================================================
namespace maniscalco
{
    template class msufsort<std::uint32_t>;
    template class msufsort<std::uint64_t>;
}
