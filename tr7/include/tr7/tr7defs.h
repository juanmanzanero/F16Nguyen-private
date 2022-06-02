#ifndef TR7_TR7DEFS_H
#define TR7_TR7DEFS_H
#pragma once


#include <stdexcept>
#include <iostream>


//
// Defines general macros, mainly to ensure both Visual Studio and gcc
// compatibility.
//


#ifdef _MSC_VER
#define constexpr_variable const
#define constexpr_ inline
#define pure_function __forceinline
#define if_constexpr if
#define enum_class enum
#define unused_variable
#define static_inline_variable static
#define std_filesystem std::experimental::filesystem

#ifndef _DEBUG
#define NDEBUG
#endif

#else
#define constexpr_variable constexpr
#define constexpr_ constexpr
#if defined(__GNUC__) || defined(__GNUG__)
#define pure_function __attribute__((always_inline))
#else
#define pure_function
#endif
#define if_constexpr if constexpr
#define enum_class enum class
#define unused_variable [[maybe_unused]]
#define static_inline_variable static inline
#define std_filesystem std::filesystem

#ifndef NDEBUG
#define _DEBUG
#if defined(__GNUC__) || defined(__GNUG__)
#define _GLIBCXX_DEBUG
#define _GLIBCXX_DEBUG_PEDANTIC
#endif
#endif

#endif


// in true c++20, functions like std::fill, and std::transform and
// std::accumulate will be constexpr
#define future_constexpr_ inline


// a handy #define to raise exceptions
#define TR7_THROW_RUNTIME_ERROR(msg) { std::runtime_error excp(msg); \
                                       std::cerr << excp.what() << std::endl; \
                                        throw excp; }


// handy findable print
#define TR7_DEBUG_COUT(MSG) { std::cout << #MSG << " = " << MSG << std::endl; }


//
// I'll also place here some standard STL functions that
// our compilers don't have...
//

#ifdef _MSC_VER
    // VS2017 does not have std::clamp...
#include <functional>
namespace std
{
    template<typename T, class Compare = std::less<>>
    constexpr const T& clamp(const T &val, const T &lo, const T &hi, Compare comp = {})
    {
        return comp(val, lo) ? lo : comp(hi, val) ? hi : val;
    }


}

#endif


namespace std
{
    // neither gcc-8 nor VS2017 have the scan routines from <numeric>...
    template<typename InputIt, typename OutputIt, typename T, typename BinaryOperation>
    constexpr OutputIt exclusive_scan(InputIt first, InputIt last,
                                      OutputIt d_first, T init, BinaryOperation binary_op)
    {
        //
        // Exclusive-scan algorithm, see std::exclusive_scan.
        //
        // TODO: change to std::exclusive_scan (not implemented in gcc-8,
        // this is a copy from github/gcc-mirror).
        //

        while (first != last) {
            auto v = init;
            init = binary_op(init, *first);
            ++first;
            *d_first++ = std::move(v);

        }

        return d_first;

    }

    template<typename InputIt, typename OutputIt, typename BinaryOperation>
    constexpr OutputIt inclusive_scan(InputIt first, InputIt last,
                                      OutputIt d_first, BinaryOperation binary_op)
    {
        //
        // Inclusive-scan algorithm, see std::inclusive_scan.
        //
        // TODO: change to std::inclusive_scan (not implemented in gcc-8,
        // this is a copy from github/gcc-mirror).
        //

        const auto inclusive_scan_ = [](auto first, auto last,
                                        auto result, auto binary_op, auto init)
        {
            for (; first != last; ++first) {
                *result++ = init = binary_op(init, *first);
            }

            return result;

        };

        if (first != last) {
            auto init = *first;
            *d_first++ = init;
            ++first;
            if (first != last) {
                d_first = inclusive_scan_(first, last, d_first,
                                          binary_op, std::move(init));

            }
        }

        return d_first;

    }

}


#endif
