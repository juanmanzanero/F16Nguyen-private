#ifndef TR7_TR7TUPLE_H
#define TR7_TR7TUPLE_H


#include <tuple>

#include "tr7/tr7defs.h"


//
// Defines handy functions and type traits
// to operate on std::tuples.
//


namespace tr7
{
    //
    // Implementation of different variants of compile-time loops,
    // that serve mainly to traverse std::tuples.
    //

#ifndef _MSC_VER
    template<typename Functor, typename IndexType, IndexType... Is>
    constexpr pure_function void index_for_impl(Functor functor, std::integer_sequence<IndexType, Is ...>)
    {
        (functor(std::integral_constant<IndexType, Is>{}), ...);
    }

    template<std::size_t N, typename Functor, typename IndexType = std::size_t>
    constexpr pure_function void index_for(Functor f)
    {
        index_for_impl(f, std::make_integer_sequence<IndexType, N>());
    }


    template<typename Tuple, typename Functor, typename IndexType, IndexType... Is>
    constexpr pure_function void tuple_for_impl(Tuple &&tuple, Functor f, std::integer_sequence<IndexType, Is...>)
    {
        (f(std::get<Is>(std::forward<Tuple>(tuple))), ...);
    }

    template<typename Tuple, typename Functor, typename IndexType = std::size_t>
    constexpr pure_function void tuple_for(Tuple &&tuple, Functor f)
    {
        tuple_for_impl(std::forward<Tuple>(tuple), f,
            std::make_integer_sequence<IndexType, std::tuple_size_v<std::decay_t<Tuple> > >());
    }


    template<typename Tuple, typename Functor, typename IndexType, IndexType... Is>
    constexpr pure_function void tuple_index_for_impl(Tuple &&tuple, Functor f, std::integer_sequence<IndexType, Is...>)
    {
        (f(std::get<Is>(std::forward<Tuple>(tuple)), std::integral_constant<IndexType, Is>{}), ...);
    }

    template<typename Tuple, typename Functor, typename IndexType = std::size_t>
    constexpr pure_function void tuple_index_for(Tuple &&tuple, Functor f)
    {
        tuple_index_for_impl(std::forward<Tuple>(tuple), f,
            std::make_integer_sequence<IndexType, std::tuple_size_v<std::decay_t<Tuple> > >());
    }


#else
    //
    // VS2017 hasn't got fold expressions... we'll implement
    // the loops using compile-time recursion...
    //

    template<std::size_t N,
             typename Functor,
             typename IndexType = std::size_t,
             IndexType I = 0u,
             typename std::enable_if_t<I >= N>* = nullptr>
    constexpr pure_function void index_for(Functor) {}

    template<std::size_t N,
             typename Functor,
             typename IndexType = std::size_t,
             IndexType I = 0u,
             typename std::enable_if_t<I < N>* = nullptr>
    constexpr pure_function void index_for(Functor f)
    {
        f(std::integral_constant<IndexType, I>{});
        index_for<N, Functor, IndexType, I + 1>(f);
    }


    template<typename Tuple, typename Functor,
             typename IndexType = std::size_t,
             IndexType I = 0u,
             typename std::enable_if_t<I >= std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
    constexpr pure_function void tuple_for(Tuple &&, Functor) {}

    template<typename Tuple, typename Functor,
             typename IndexType = std::size_t,
             IndexType I = 0u,
             typename std::enable_if_t<I < std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
    constexpr pure_function void tuple_for(Tuple &&tuple, Functor f)
    {
        f(std::get<I>(tuple));
        tuple_for<Tuple, Functor, IndexType, I + 1>(tuple, f);
    }


    template<typename Tuple, typename Functor,
             typename IndexType = std::size_t,
             IndexType I = 0u,
             typename std::enable_if_t<I >= std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
    constexpr pure_function void tuple_index_for(Tuple &&, Functor) {}

    template<typename Tuple, typename Functor,
             typename IndexType = std::size_t,
             IndexType I = 0u,
             typename std::enable_if_t<I < std::tuple_size_v<std::decay_t<Tuple> > >* = nullptr>
    constexpr pure_function void tuple_index_for(Tuple &&tuple, Functor f)
    {
        f(std::get<I>(tuple), std::integral_constant<IndexType, I>{});
        tuple_index_for<Tuple, Functor, IndexType, I + 1>(tuple, f);
    }


#endif


    //
    // Performs std::decay on every member type of an std::tuple.
    //

    template <typename ... Ts>
    constexpr auto decay_types_of_a_tuple(const std::tuple<Ts...> &) ->
        std::tuple<std::remove_cv_t<std::remove_reference_t<Ts> >... >;

    template <typename Tuple>
    using decay_tuple = decltype(decay_types_of_a_tuple(std::declval<Tuple>()));


}


#endif
