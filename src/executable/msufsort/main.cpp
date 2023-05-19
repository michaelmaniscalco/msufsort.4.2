
#include <library/msufsort.h>

#include <fmt/format.h>

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <span>
#include <thread>

#include <concepts>



namespace
{

    //=========================================================================
    std::vector<char> load_file
    (
        std::string const & path
    )
    {
     //   return {'a', 'b', 'a', 'b', 'a', 'b', 'a', 'b'};
        std::vector<char> destination;
        std::ifstream source(path, std::ios_base::in | std::ios_base::binary);
        if (!source)
            throw std::runtime_error(fmt::format("Failed to load file: {}", path));
        source.seekg(0, std::ios_base::end);
        destination.resize(source.tellg());
        source.seekg(0, std::ios_base::beg);
        source.read(destination.data(), destination.size());
        source.close();
        return destination;
    }

/*
//==============================================================================
inline void insertion_sort
(
    int * begin,
    int * end
)
{
    std::int32_t partitionSize = std::distance(begin, end);
    if (partitionSize <= 2)
    {
        if ((partitionSize == 2) && (begin[0] > begin[1]))
            std::swap(begin[0], begin[1]);
        return;
    }

    int value[8];
    auto n = 0;
    for (auto cur = begin; cur < end; ++cur)
        value[n++] = cur[0];

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
auto partition
(
    std::span<int> suffixArray
) -> std::pair<std::span<int>, std::span<int>>
{
    // less than p1 | equal p1 | between p1 and p2 | equal p2 | greater than p2
    //              a          b         c         d          e
    auto a = suffixArray.data();
    auto b = a;
    auto d = suffixArray.data() + suffixArray.size() - 1;
    auto e = d;
    auto leftPivot = *b;
    auto rightPivot = *d;
    if (leftPivot > rightPivot)
        std::swap(*b, *d), std::swap(leftPivot, rightPivot);
    ++b;
    --d;
    auto c = b;
    while (c < d)
    {
        std::swap(*c, *d);
        auto value = *c;
        while ((c <= d) && (value < rightPivot))
        {
            if (value <= leftPivot)
            {
                std::swap(*b, *c);
                if (value < leftPivot)
                    std::swap(*b, *a++);
                ++b;
            }
            value = *++c;
        }

        value = *d;
        while ((c <= d) && (value >= rightPivot))
        {
            if (value > rightPivot)
                std::swap(*d, *e--);
            value = *--d;
        }
    }
    ++d;
    ++e;
    return {{a, std::distance(a, b)}, {d, std::distance(d, e)}};
}


void quick_sort
(
    std::span<int> input
)
{
    std::vector<std::span<int>> stack;
    stack.push_back(input);
    while (!stack.empty())
    {
        auto data = stack.back();
        stack.pop_back();
        if (data.size() < 4)
            insertion_sort(data.data(), data.data() + data.size());
        else
        {
            auto [p1, p2] = partition(data);
            auto size = (p1.begin() - data.begin());
            if (size)
                stack.push_back({data.begin(), size});
        //    if (p1.size())
        //        stack.push_back(p1);
            size = (p2.begin() - p1.end());
            if (size)
                stack.push_back({p1.end(), size});
        //    if (p2.size())
         //       stack.push_back(p2);
            size = (data.end() - p2.end());
            if (size)
                stack.push_back({p2.end(), size});
        }
    }
}
*/
} // namespace


//=============================================================================
int main
(
    int argc,
    char ** argv
)
{
    try
    {
        using namespace maniscalco;
        std::uint32_t threadCount{1};
        if (argc > 2)
        {
            try
            {
                threadCount = std::uint32_t(std::atoi(argv[2]));
            }
            catch (...)
            {
                std::cerr << "invalid thread count: " << argv[2] << "\n";
            }
        }
        if (threadCount > std::thread::hardware_concurrency())
        {
            threadCount = std::uint32_t(std::thread::hardware_concurrency());
            std::cerr << "restricting thread count to max concurrency of " << threadCount << '\n';
        }
        std::cout << "concurrency = " << threadCount << "\n";
        auto input = load_file(argv[1]);
        make_suffix_array<decltype(input), std::uint32_t>(input, threadCount);
    }
    catch (std::exception const & exception)
    {
        std::cout << exception.what() << '\n';
    }

    return 0;
}
