#ifndef TR7_TR7AD_H
#define TR7_TR7AD_H
#pragma once


#include <tuple>
#include <functional>


#ifdef TR7_WITH_AD
#include "cppad/cppad.hpp"
#endif


//
// Contains all of the shamefule tricks we need
// to make our codes work with autodiff. We only
// support COIN-OR's "CppAD" library.
//


namespace tr7
{
    //
    // a type_trait to identify if a scalar number
    // is of AD-type
    //

    template<typename T>
    struct is_AD_scalar : std::false_type {};

#ifdef TR7_WITH_AD
    template<typename T>
    struct is_AD_scalar<CppAD::AD<T> > : std::true_type {};

#endif

    template<typename T>
    inline constexpr bool is_AD_scalar_v = is_AD_scalar<T>::value;


    //
    // std::true_type if ANY of the member types
    // of an std::tuple is an AD scalar
    //

    template<typename Tuple, typename = void>
    struct any_tuple_type_is_AD_scalar : std::false_type {};

    template<typename... T>
    struct any_tuple_type_is_AD_scalar<std::tuple<T... >,
        std::void_t<std::enable_if_t<std::disjunction_v<is_AD_scalar<typename std::decay_t<T> >... > > > > : std::true_type {};


    //
    // std::true_type if ALL of the member types
    // of an std::tuple are AD scalars
    //

    template<typename Tuple, typename = void>
    struct all_tuple_types_are_AD_scalars : std::false_type {};

    template<typename... T>
    struct all_tuple_types_are_AD_scalars<std::tuple<T... >,
        std::void_t<std::enable_if_t<std::conjunction_v<is_AD_scalar<typename std::decay_t<T> >... > > > > : std::true_type {};


    //
    // Force a scalar value of an AD-type to become non-AD
    // (you SHOULD know the consequences...)
    //

    template<typename T>
    constexpr pure_function auto AD_force_scalar_value(const T &value)
    {
#ifdef TR7_WITH_AD
        if constexpr (!is_AD_scalar_v<T>) {
#endif
            return value;
#ifdef TR7_WITH_AD
        }
        else {
            return CppAD::Value(CppAD::Var2Par(value));
        }
#endif
    }


//
// "constexpr" variables are not compatible
// with AD...
//

#ifdef TR7_WITH_AD
#ifdef constexpr_variable
#undef constexpr_variable
#define constexpr_variable const
#endif
#endif


//
// Declare <cmath> operations so that we have
// them well located.
//

#ifdef TR7_WITH_AD
#define DECLARE_TR7_AD_CMATHFUN(CMATHFUN) template<typename... Args> \
                                          constexpr pure_function auto AD_##CMATHFUN(Args&&... args) \
                                          { \
                                              if constexpr (!any_tuple_type_is_AD_scalar<std::tuple<Args...> >::value) { \
                                                  return std::CMATHFUN(std::forward<Args>(args)...); \
                                              } \
                                              else { \
                                                  return CMATHFUN(args...); \
                                              } \
                                          }

#else
#define DECLARE_TR7_AD_CMATHFUN(CMATHFUN) template<typename... Args> \
                                          constexpr pure_function auto AD_##CMATHFUN(Args&&... args) \
                                          { \
                                              return std::CMATHFUN(std::forward<Args>(args)...); \
                                          }

#endif

DECLARE_TR7_AD_CMATHFUN(abs)
DECLARE_TR7_AD_CMATHFUN(sqrt)
DECLARE_TR7_AD_CMATHFUN(pow)
DECLARE_TR7_AD_CMATHFUN(exp)
DECLARE_TR7_AD_CMATHFUN(log)
DECLARE_TR7_AD_CMATHFUN(log10)
DECLARE_TR7_AD_CMATHFUN(sin)
DECLARE_TR7_AD_CMATHFUN(cos)
DECLARE_TR7_AD_CMATHFUN(tan)
DECLARE_TR7_AD_CMATHFUN(asin)
DECLARE_TR7_AD_CMATHFUN(acos)
DECLARE_TR7_AD_CMATHFUN(atan)
DECLARE_TR7_AD_CMATHFUN(atan2)
#undef DECLARE_TR7_AD_CMATHFUN


    //
    // An equivalent of std::clamp that is AD-friendly.
    //

    template<typename T, class Compare = std::less<>>
    constexpr pure_function const T& AD_clamp(const T &val, const T &lo, const T &hi, Compare comp = {})
    {
#ifdef TR7_WITH_AD
        if constexpr (!is_AD_scalar_v<T>) {
#endif
            return std::clamp(val, lo, hi, comp);
#ifdef TR7_WITH_AD
        }
        else {
            return comp(val, lo) ? lo : comp(hi, val) ? hi : val;
        }
#endif
    }


}


#endif