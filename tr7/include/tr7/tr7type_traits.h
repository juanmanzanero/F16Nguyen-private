#ifndef TR7_TR7TYPE_TRAITS_H
#define TR7_TR7TYPE_TRAITS_H
#pragma once


#include <type_traits>
#include <iterator>
#include <complex>
#include <array>


//
// Defines some useful type_traits.
//


namespace tr7
{
    //
    // std::true_type if the tested type is an
    // specialization of a template class.
    //

    template<typename T, template<typename...> class>
    struct is_specialization : std::false_type {};

    template<template<typename...> class T, typename... Ts>
    struct is_specialization<T<Ts...>, T> : std::true_type {};


    //
    // std::true_type if the tested type is floating
    // point or complex.
    //

    template<typename T>
    struct is_floating_point_or_complex : std::is_floating_point<T> {};
    template<typename T>
    struct is_floating_point_or_complex<std::complex<T>> : std::is_floating_point<T> {};


    //
    // std::true_type if the tested type has a ".clear()"
    // member function.
    //

    template<typename T, typename = void>
    struct has_clear_member_fun : std::false_type {};

    template<typename T>
    struct has_clear_member_fun<T,
                                std::void_t<decltype(std::declval<T &>().clear())> > : std::true_type {};


    //
    // std::true_type if the tested type has a ".resize()"
    // member function.
    //

    template<typename T, typename = void>
    struct has_resize_member_fun : std::false_type {};

    template<typename T>
    struct has_resize_member_fun<T,
                                 std::void_t<decltype(std::declval<T &>().resize(std::size_t{}))> > : std::true_type {};


    //
    // std::true_type if the tested type has a ".size()"
    // member function.
    //

    template<typename T, typename = void>
    struct has_size_member_fun : std::false_type {};

    template<typename T>
    struct has_size_member_fun<T,
                               std::void_t<decltype(std::declval<T &>().size())> > : std::true_type {};


    //
    // A helper to find the type that we get
    // when applying the bracket operator
    // to an object (e.g., "decltype(v[i])").
    //

    template<class C, typename IndexType = std::size_t>
    using typeof_bracket_operator = std::decay_t<decltype(std::declval<C &>()[std::declval<IndexType>()])>;


    //
    // std::true_type if the tested type represents
    // an iterator/pointer.
    //

    template<typename It, typename = void>
    struct is_iterator : std::false_type {};

    template<typename It>
    struct is_iterator<It,
                       std::void_t<std::enable_if_t<!std::is_same_v<typename std::iterator_traits<It>::value_type, void> > > > : std::true_type {};


    //
    // std::true_type if the tested type represents
    // a const iterator/pointer.
    //

    template<typename It, typename = void>
    struct is_const_iterator : std::false_type {};

    template<typename It>
    struct is_const_iterator<It,
                             std::void_t<std::enable_if_t<std::is_const_v<std::remove_pointer_t<typename std::iterator_traits<It>::pointer> > > > > : std::true_type {};


    //
    // std::true_type if the tested type is
    // an std::array.
    //

    template<typename Arr, typename = void>
    struct is_std_array : std::false_type {};

    template<typename T, std::size_t N>
    struct is_std_array<std::array<T, N> > : std::true_type {};


    //
    // std::common_type applied to the member
    // types of an std::tuple
    //

    template<typename Tuple>
    struct tuple_common_type {};

    template<typename... T>
    struct tuple_common_type<std::tuple<T... > > : std::common_type<T...> {};

}


#endif
