#ifndef F16_NGUYEN_CONTIGUOUS_RANGES
#define F16_NGUYEN_CONTIGUOUS_RANGES
#pragma once


#include "tr7/tr7defs.h"
#include "tr7/tr7type_traits.h"


//
// Defines the "*contiguous_*_range" family
// of classes, which help us to check sizes
// and bounds when working with pointers/containers
// to which we apply the bracket operator (i.e., "v[i]").
//


namespace F16_Nguyen
{
    //
    // "contiguous_range" and "const_contiguous_range",
    // to define a range whose size is not known at
    // compile time.
    //

    template<typename It>
    class contiguous_range : public std::iterator_traits<It>
    {
        //
        // Represents a contiguous range of elements
        // that are defined by a "begin" iterator/pointer
        // and a size, such that the elements in
        // "[begin, begin + size)" can be accessed
        // through the bracket operator.
        //

    public:

        using base_type = std::iterator_traits<It>;
        using typename base_type::reference;

        using iterator = It;


        constexpr contiguous_range(It begin, std::size_t size) : _begin{ begin }, _size{ size } {}

        constexpr contiguous_range(It begin, It end) : contiguous_range(begin, std::distance(begin, end)) {}

        constexpr iterator begin() const { return _begin; }
        constexpr iterator end() const { return _begin + _size; }

        constexpr const std::size_t& size() const { return _size; }

        template<typename SizeType>
        constexpr reference operator[](const SizeType &pos) const { return *std::next(_begin, pos); }

        template<typename SizeType>
        constexpr reference at(const SizeType &pos) const noexcept(false)
        {
            if (static_cast<std::size_t>(pos) < _size) {
                return this->operator[](pos);
            }
            else {
                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::contiguous_range::at: out of bounds.")
            }

        }


    private:

        It _begin;
        std::size_t _size;

    };


    namespace detail
    {
        //
        // A helper type trait to convert non-const iterator
        // types into pointer-to-const types.
        //

        template<typename It, typename = void>
        struct force_const_iterator {};

        template<typename It>
        struct force_const_iterator<It,
                                   std::void_t<std::enable_if_t<tr7::is_iterator<It>::value &&
                                                                !tr7::is_const_iterator<It>::value > > >
        {
            using type = const typename std::iterator_traits<It>::value_type *;
        };

        template<typename It>
        struct force_const_iterator<It,
                                  std::void_t<std::enable_if_t<tr7::is_iterator<It>::value &&
                                                               tr7::is_const_iterator<It>::value > > >
        {
            using type = It;
        };

    }

    template<typename It>
    struct const_contiguous_range : contiguous_range<typename detail::force_const_iterator<It>::type>
    {
        //
        // The "const iterator" version of "contiguous_range",
        // i.e., a range in which we can only access const
        // elements. If template parameter "It" points to
        // non-const values, we'll convert it into a
        // pointer-to-const type.
        //

        using base_type = contiguous_range<typename detail::force_const_iterator<It>::type>;


        template<typename It_,
                 typename... Args,
                 typename std::enable_if_t<!tr7::is_const_iterator<It_>::value>* = nullptr>
        constexpr const_contiguous_range(It_ begin, Args&&... args) : base_type(&(*begin), std::forward<Args>(args)...) {}

        template<typename It_,
                 typename... Args,
                 typename std::enable_if_t<tr7::is_const_iterator<It_>::value>* = nullptr>
        constexpr const_contiguous_range(It_ begin, Args&&... args) : base_type(begin, std::forward<Args>(args)...) {}

    };


    //
    // "static_contiguous_range" and
    // "const_static_contiguous_range",
    // to define a range whose size is known
    // at compile time.
    //

    template<std::size_t Size, typename It>
    class static_contiguous_range : public std::iterator_traits<It>
    {
        //
        // Defines a range like "F16_Nguyen::contiguous_range",
        // but in this case the size is a compile-time constant.
        //

    public:

        using base_type = std::iterator_traits<It>;
        using typename base_type::reference;

        using iterator = It;

        static constexpr std::size_t _size = Size;


        constexpr static_contiguous_range(It begin) : _begin{ begin } {}

        constexpr iterator begin() const { return _begin; }
        constexpr iterator end() const { return _begin + _size; }

        static constexpr const std::size_t& size() { return _size; }


        template<typename SizeType>
        constexpr reference operator[](const SizeType &pos) const { return *std::next(_begin, pos); }

        template<typename SizeType>
        constexpr reference at(const SizeType &pos) const noexcept(false)
        {
            if (static_cast<std::size_t>(pos) < _size) {
                return this->operator[](pos);
            }
            else {
                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::static_contiguous_range::at: out of bounds.")
            }

        }

        template<std::size_t Pos>
        constexpr reference get() const
        {
            static_assert(Pos < _size, "F16_Nguyen::static_contiguous_range::get: out of bounds.");
            return this->operator[](Pos);
        }


    private:

        It _begin;

    };


    template<std::size_t Size,
             typename It>
    struct const_static_contiguous_range : static_contiguous_range<Size, typename detail::force_const_iterator<It>::type>
    {
        //
        // The "const iterator" version of "static_contiguous_range",
        // i.e., a range in which we can only access const
        // elements.
        //

        using base_type = static_contiguous_range<Size, typename detail::force_const_iterator<It>::type>;


        template<typename It_,
                 typename std::enable_if_t<!tr7::is_const_iterator<It_>::value>* = nullptr>
        constexpr const_static_contiguous_range(It_ begin) : base_type(&(*begin)) {}

        template<typename It_,
                 typename std::enable_if_t<tr7::is_const_iterator<It_>::value>* = nullptr>
        constexpr const_static_contiguous_range(It_ begin) : base_type(begin) {}

    };


    template<class RangeType, typename = void>
    struct is_static_contiguous_range : std::false_type {};

    template<std::size_t Size, typename It>
    struct is_static_contiguous_range<static_contiguous_range<Size, It> > : std::true_type {};

    template<std::size_t Size, typename It>
    struct is_static_contiguous_range<const_static_contiguous_range<Size, It> > : std::true_type {};

    template<class RangeType, typename = void>
    struct is_const_static_contiguous_range : std::false_type {};

    template<std::size_t Size, typename It>
    struct is_const_static_contiguous_range<const_static_contiguous_range<Size, It> > : std::true_type {};


    //
    // The "makers" of ranges, use these ones instead of
    // the ctors so that the "It" template parameter can
    // get auto-deduced.
    //

    template<typename It, typename SizeOrEndType>
    constexpr auto make_contiguous_range(It begin, SizeOrEndType size_or_end)
    {
        return contiguous_range<It>(begin, size_or_end);
    }

    template<typename It, typename SizeOrEndType>
    constexpr auto make_const_contiguous_range(It begin, SizeOrEndType size_or_end)
    {
        return const_contiguous_range<It>(begin, size_or_end);
    }

    template<std::size_t Size,
             typename It>
    constexpr auto make_static_contiguous_range(It begin)
    {
        return static_contiguous_range<Size, It>(begin);
    }

    template<std::size_t Size,
             typename It>
    constexpr auto make_const_static_contiguous_range(It begin)
    {
        return const_static_contiguous_range<Size, It>(begin);
    }


    //
    // A collection of functions to check sizes of ranges. They
    // do nothing if their size cannot be calculated.
    //

    template<class RangeType>
    constexpr void ensure_size(std::size_t required_size, const RangeType &r) noexcept(false)
    {
        if constexpr (tr7::has_size_member_fun<RangeType>::value) {
            //assert(r.size() == required_size)
            if (r.size() != required_size) {
                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::ensure_size: the input range doesn't have the required size.")
            }

        }

    }

    template<std::size_t RequiredSize,
             class RangeType>
    constexpr void ensure_size(const RangeType &r)
    {
        if constexpr (is_static_contiguous_range<RangeType>::value) {
            static_assert(RangeType::_size == RequiredSize,
                          "F16_Nguyen::ensure_size: the input range doesn't have the required size.");

        }
        else if constexpr (tr7::is_std_array<RangeType>::value ||
                           tr7::is_specialization<RangeType, std::tuple>::value ||
                           tr7::is_specialization<RangeType, std::pair>::value) {

            static_assert(std::tuple_size_v<RangeType> == RequiredSize,
                          "F16_Nguyen::ensure_size: the input range doesn't have the required size.");

        }
        else {
            ensure_size(RequiredSize, r);

        }

    }

    template<class RangeType>
    constexpr void ensure_minimum_size(std::size_t required_minimum_size, const RangeType &r) noexcept(false)
    {
        if constexpr (tr7::has_size_member_fun<RangeType>::value) {
            //assert(r.size() >= required_minimum_size)
            if (r.size() < required_minimum_size) {
                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::ensure_minimum_size: the input range doesn't reach the required minimum size.")
            }

        }

    }

    template<std::size_t RequiredMinimumSize,
             class RangeType>
    constexpr void ensure_minimum_size(const RangeType &r)
    {
        if constexpr (is_static_contiguous_range<RangeType>::value) {
            static_assert(RangeType::_size >= RequiredMinimumSize,
                          "F16_Nguyen::ensure_minimum_size: the input range doesn't reach the required minimum size.");

        }
        else if constexpr (tr7::is_std_array<RangeType>::value ||
                           tr7::is_specialization<RangeType, std::tuple>::value ||
                           tr7::is_specialization<RangeType, std::pair>::value) {

            static_assert(std::tuple_size_v<RangeType> >= RequiredMinimumSize,
                          "F16_Nguyen::ensure_minimum_size: the input range doesn't reach the required minimum size.");

        }
        else {
            ensure_minimum_size(RequiredMinimumSize, r);

        }

    }


}


#endif
