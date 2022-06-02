#ifndef TR7_LOOKUP_H
#define TR7_LOOKUP_H
#pragma once


#include <iostream>
#include <type_traits>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>

#include "tr7/tr7type_traits.h"
#include "tr7/tr7tuple.h"


//
// Defines the "Prelookup" and "Lookup_table" classes, to interpolate on
// a sorted grid.
//


namespace tr7
{
    struct lookup
    {
        //
        // Holds some options & setup to perform lookup operations.
        //

        // enumeration of available interpolation methods. NOTE: the "spline"
        // method is plagiarized from Simulink's "Lookup Table" block, which is
        // MUCH SLOWER than Matlab's built-in "griddedInterpolant" class
        enum class interpolation{ linear, nearest, previous, next, spline };


        // enumeration of available extrapolation methods: with "nearest", we
        // take the last values in the table, with "allow", we extrapolate using
        // the selected interpolation method, with "linear" we force linear
        // extrapolation. NOTE: if "interpolation = nearest/previous/next", then
        // "extrapolation = linear" is NOT SUPPORTED, and the two remaining methods
        // have exactly the same behaviour
        enum class extrapolation{ allow, nearest, linear };


        // search (a.k.a. "prelookup") methods: "linear search", "binary search"
        // (i.e., bisection), "hunt-locate" (a combination of the previous two),
        // "equispaced" (only valid if all of the grid vectors have constant spacings)
        enum class prelookup{ linsearch, binsearch, hunt_locate, equispaced };


        // the type we'll use as index in our lookup operations
        using index_type = std::size_t;


        struct extrapolation_flag
        {
            //
            // Enumerates the directions in which an extrapolation
            // can take place.
            //
            // NOTE: DON'T TOUCH THESE NUMERICAL VALUES!!
            //

            static constexpr index_type lower = 0u;
            static constexpr index_type upper = 1u;
            static constexpr index_type none = 2u;

        };

    };


    template<lookup::prelookup PrelookupOption = lookup::prelookup::binsearch,
             bool UsePreviousIndex = false,
             typename GridVectorType = std::vector<double> >
    class Prelookup
    {
        //
        // Represents the "prelookup" (i.e., index search) operation
        // that we have to carry out in all of a lookup table's
        // grid vectors.
        //

    public:

        using index_type = typename lookup::index_type;
        using signed_index_type = std::ptrdiff_t;
        using grid_vector_type = GridVectorType;
        using prelookup_result = std::array<index_type, 2>;

        static constexpr lookup::prelookup prelookup_option = PrelookupOption;
        static constexpr bool use_previous_index = UsePreviousIndex;


        constexpr Prelookup() = default;

        template<typename GridVector>
        Prelookup(GridVector &&grid_vector_)
        {
            grid_vector(std::forward<GridVector>(grid_vector_));
        }


        constexpr const grid_vector_type& grid_vector() const { return _grid_vector; }

        template<typename GridVector>
        inline void grid_vector(GridVector &&grid_vector)
        {
            grid_vector_impl(std::forward<GridVector>(grid_vector));
        }


        constexpr const index_type& max_index() const { return _max_index; }


        constexpr bool empty() const
        {
            return _grid_vector.empty();
        }

        constexpr void clear()
        {
            _max_index = 0u;
            _extra_data = decltype(_extra_data){};

#ifndef _MSC_VER
            if constexpr (has_clear_member_fun<grid_vector_type>::value) {
#else
            if constexpr (is_specialization<grid_vector_type, std::vector>::value ||
                          is_specialization<grid_vector_type, std::map>::value ||
                          is_specialization<grid_vector_type, std::unordered_map>::value) {
#endif
                _grid_vector.clear();

            }
            else {
                _grid_vector = grid_vector_type{};

            }

        }


        template<typename QueryType>
        constexpr prelookup_result operator()(const QueryType &q) const
        {
            //
            // Prelookup operation, protected for extrapolation.
            //

            prelookup_result ret{};
            auto &idx_lo = ret[0];
            auto &e = ret[1];
            if (q > _grid_vector[_max_index]) {
                idx_lo = _max_index - index_type{ 1u };
                e = lookup::extrapolation_flag::upper;

                if constexpr (prelookup_option == lookup::prelookup::hunt_locate) {
                    _extra_data.cor = false;
                }

            }
            else if (q < _grid_vector[0]) {
                idx_lo = index_type{ 0u };
                e = lookup::extrapolation_flag::lower;

                if constexpr (prelookup_option == lookup::prelookup::hunt_locate) {
                    _extra_data.cor = false;
                }

            }
            else {
                idx_lo = prelookup_inside(q);
                e = lookup::extrapolation_flag::none;

            }

            if constexpr (use_previous_index) {
                _extra_data.prev_index = idx_lo;
            }

            return ret;

        }


        template<bool RequiresLinearFraction,
                 lookup::extrapolation ExtrapolationOption = lookup::extrapolation::allow,
                 typename QueryType>
        constexpr auto evaluate(const QueryType &q) const
        {
            //
            // Evaluates the prelookup operation and optionally returns
            // the associated linear fraction (if template parameter
            // "CalculateLinearFraction" is true), whose formula
            // depends on the extrapolation method we're using (i.e.,
            // on the value of the "ExtrapolationOption" parameter).
            //

            if constexpr (RequiresLinearFraction) {
                const auto p = this->operator()(q);
                return std::make_pair(p, linear_frac<ExtrapolationOption>(q, p));

            }
            else {
                return this->operator()(q);

            }

        }


    private:

        grid_vector_type _grid_vector;
        index_type _max_index;


        template<bool, typename = void>
        struct basic_extra_data_prelookup
        {
            template<typename... Args>
            static constexpr void build(Args&&...) {}
        };

        template<bool UsePreviousIndex_>
        struct basic_extra_data_prelookup<UsePreviousIndex_,
                                          std::void_t<std::enable_if_t<UsePreviousIndex_> > >
        {
            template<typename... Args>
            constexpr void build(index_type max_index, Args&&...)
            {
                // initialize the prev_indices to the midpoint of each grid_vector, for example
                // NOTE: x >> 1u == floor(x / 2)
                prev_index = max_index >> 1u;

            }


            index_type prev_index{ 0u };

        };

        template<lookup::prelookup PrelookupOption_,
                 bool UsePreviousIndex_, typename = void>
        struct extra_data_prelookup : basic_extra_data_prelookup<UsePreviousIndex_> {};

        template<lookup::prelookup PrelookupOption_, bool UsePreviousIndex_>
        struct extra_data_prelookup <PrelookupOption_, UsePreviousIndex_,
                                     std::void_t<std::enable_if_t<PrelookupOption_ == lookup::prelookup::hunt_locate> > > :
               basic_extra_data_prelookup<UsePreviousIndex_>
        {
            template<typename... Args>
            constexpr void build(index_type max_index, Args&&...)
            {
                basic_extra_data_prelookup<UsePreviousIndex_>::build(max_index);

                // initialize the "cor" to false, and the "cor_dist" to
                // the typical distance for the grid vector (grid_vector_size^0.25)
                cor = false;
                cor_dist = static_cast<signed_index_type>(std::max(1.,
                                                                   std::pow(static_cast<double>(max_index) + 1., 0.25)));

            }


            bool cor{ false };
            signed_index_type cor_dist{ 0 };

        };

        template<lookup::prelookup PrelookupOption_, bool UsePreviousIndex_>
        struct extra_data_prelookup<PrelookupOption_, UsePreviousIndex_,
                                    std::void_t<std::enable_if_t<PrelookupOption_ == lookup::prelookup::equispaced> > >
        {
            using grid_type = typename grid_vector_type::value_type;

            static_assert(!UsePreviousIndex_, "tr7::Prelookup::extra_data_prelookup: the \"equispaced\" prelookup method"
                                              " can NEVER use the previous index.");


            template<typename GridVector>
            inline void build(index_type max_index, const GridVector &grid_vector) noexcept(false)
            {
                auto gv0 = std::cbegin(grid_vector);
                auto gv1 = std::next(gv0);
                const auto first_spacing = *gv1 - *gv0;

                const auto gvlast = std::prev(std::cend(grid_vector));
                const auto tol_spacing = (*gvlast - *gv0) * std::numeric_limits<grid_type>::epsilon();

                // ensure that the grid vector is actually equispaced, we only
                // admit an "extent_grid_vector * eps" discrepancy. At the end
                // we'll take the mean spacing with the intent of minimizing
                // the roundoff error
                spacing = first_spacing;
                while (gv1 != gvlast) {
                    const auto s = *++gv1 - *++gv0;
                    if (std::abs(s - first_spacing) > tol_spacing) {
                        TR7_THROW_RUNTIME_ERROR("tr7::Prelookup::extra_data_prelookup::build: input grid"
                                                " vector must be equispaced.");

                    }
                    else {
                        spacing += s;

                    }

                }

                spacing /= static_cast<grid_type>(std::distance(std::cbegin(grid_vector), gvlast));

            }


            grid_type spacing{ 0 };

        };

        mutable extra_data_prelookup<prelookup_option, use_previous_index> _extra_data;


        template<typename GridVector>
        inline void grid_vector_impl(GridVector &&grid_vector) noexcept(false)
        {
            //
            // Validates the grid and generates all
            // of the necessary data.
            //

            const auto grid_vector_size = grid_vector.size();

            // support for the empty object
            if (grid_vector_size == 0u) {
                clear();
                return;

            }

            // check that the grid vector has an adequate size
            if (grid_vector_size < 2u) {
                TR7_THROW_RUNTIME_ERROR("tr7::Prelookup::build_grid: interpolation requires"
                                        " at least two sample points in each dimension.");

            }

            // check that grid vectors are strictly monotonically increasing
            if (!std::is_sorted(std::cbegin(grid_vector), std::cend(grid_vector),
                                [](auto g1, auto g0) { return g1 <= g0; })) {

                TR7_THROW_RUNTIME_ERROR("tr7::Prelookup::build_grid: grid vectors must"
                                        " always be strictly monotonically increasing.");

            }

            // calculate max index and the extra data that depends on
            // the method we'll employ, and assign the grid vector
            _max_index = static_cast<index_type>(grid_vector_size - 1u);

            _extra_data.build(_max_index, grid_vector);

            if constexpr (std::is_same_v<GridVector, grid_vector_type>) {
                _grid_vector = std::forward<GridVector>(grid_vector);
            }
            else {
                _grid_vector = grid_vector_type(std::cbegin(grid_vector), std::cend(grid_vector));
            }

        }


        template<lookup::prelookup PrelookupOption_ = prelookup_option,
                 typename QueryType,
                 typename std::enable_if_t<PrelookupOption_ == lookup::prelookup::linsearch>* = nullptr>
        constexpr index_type prelookup_inside(const QueryType &q) const
        {
            //
            // Prelookup inside the grid with the linear search method.
            //

            index_type idx_lo;
            if constexpr (use_previous_index) {
                idx_lo = _extra_data.prev_index;
            }
            else {
                idx_lo = _max_index >> 1u;
            }

            for (; q < _grid_vector[idx_lo] && idx_lo != index_type{ 0u }; idx_lo--) {}

            const auto mm1 = _max_index - index_type{ 1u };
            while (q >= _grid_vector[idx_lo + 1] && idx_lo != mm1) {
                idx_lo++;
            }

            return idx_lo;

        }

        template<lookup::prelookup PrelookupOption_ = prelookup_option,
                 typename QueryType,
                 typename std::enable_if_t<PrelookupOption_ == lookup::prelookup::binsearch>* = nullptr>
        constexpr index_type prelookup_inside(const QueryType &q) const
        {
            //
            // Prelookup inside the grid with the binary search method.
            //

            // complete algorithm, based on "std::lower_bound"
            // (see https://en.cppreference.com/w/cpp/algorithm/lower_bound),
            // it does not admit a starting index...
            //
            //auto count = static_cast<signed_index_type>(_max_index);
            //auto idx_lo = count - signed_index_type{ 1 };
            //
            //while (count > 0) {
            //    const auto step = count >> 1u;
            //    const auto i = idx_lo - step;
            //    if (_grid_vector[i] > q) {
            //        idx_lo = i - signed_index_type{ 1 };
            //        count -= (step + signed_index_type{ 1 });
            //    }
            //    else {
            //        count = step;
            //    }
            //}
            //return idx_lo > signed_index_type{ 0 } ? static_cast<index_type>(idx_lo) : 0u;


            // the algorithm that Simulink uses
            if constexpr (use_previous_index) {
                // NOTE: this one only works if it has been previously
                // protected from extrapolation!! We also have to add
                // the extra possibility of "q == _grid_vector[_max_index]"
                if (q == _grid_vector[_max_index]) {
                    return _max_index - index_type{ 1u };
                }

                auto idx_lo = _extra_data.prev_index;
                auto i_left = index_type{ 0u };
                auto i_right = _max_index;
                bool found = false;
                while (!found) {
                    if (q < _grid_vector[idx_lo]) {
                        i_right = idx_lo - index_type{ 1u };
                        idx_lo = (i_right + i_left) >> 1u;

                    }
                    else if (q < _grid_vector[idx_lo + 1u]) {
                        found = true;

                    }
                    else {
                        i_left = idx_lo + index_type{ 1u };
                        idx_lo = (i_right + i_left) >> 1u;

                    }

                }

                return idx_lo;

            }
            else {
                index_type i_mid = _max_index >> 1u;
                auto idx_lo = index_type{ 0u };
                index_type i_right = _max_index;
                while (i_right - idx_lo > index_type{ 1u }) {
                    if (q < _grid_vector[i_mid]) {
                        i_right = i_mid;
                    }
                    else {
                        idx_lo = i_mid;
                    }

                    i_mid = (i_right + idx_lo) >> 1u;

                }

                return idx_lo;

            }

        }

        template<lookup::prelookup PrelookupOption_ = prelookup_option,
                 typename QueryType,
                 typename std::enable_if_t<PrelookupOption_ == lookup::prelookup::hunt_locate>* = nullptr>
        constexpr index_type prelookup_inside(const QueryType &q) const
        {
            //
            // Prelookup inside the grid with the "hunt-locate" method.
            //

            // hunt phase
            index_type starting_index;
            if constexpr (use_previous_index) {
                starting_index = _extra_data.prev_index;
            }
            else {
                starting_index = _max_index >> 1u;
            }

            auto idx_lo = starting_index;
            index_type idx_up;
            if (_extra_data.cor) {
                auto inc = index_type{ 1u };

                if (q >= _grid_vector[idx_lo]) {
                    // hunt up
                    while (true) {
                        idx_up = idx_lo + inc;
                        if (idx_up >= _max_index) {
                            idx_up = _max_index;
                            break;

                        }
                        else if (q < _grid_vector[idx_up]) {
                            break;

                        }
                        else {
                            idx_lo = idx_up;
                            inc += inc;

                        }

                    }

                }
                else {
                    // hunt down
                    idx_up = idx_lo;

                    while (true) {
                        if (inc >= idx_up) {
                            idx_lo = index_type{ 0u };
                            break;

                        }
                        else {
                            idx_lo = idx_up - inc;

                        }

                        if (q >= _grid_vector[idx_lo]) {
                            break;

                        }
                        else {
                            idx_up = idx_lo;
                            inc += inc;

                        }

                    }

                }

            }
            else {
                idx_lo = index_type{ 0u };
                idx_up = _max_index;

            }

            // bisection phase
            while (idx_up - idx_lo > index_type{ 1u }) {
                const auto idx_m = (idx_lo + idx_up) >> 1u;
                if (q >= _grid_vector[idx_m]) {
                    idx_lo = idx_m;
                }
                else {
                    idx_up = idx_m;
                }

            }

            // save correlation
            _extra_data.cor = std::abs(static_cast<signed_index_type>(idx_lo) -
                                       static_cast<signed_index_type>(starting_index)) <= _extra_data.cor_dist;

            return idx_lo;

        }

        template<lookup::prelookup PrelookupOption_ = prelookup_option,
                 typename QueryType,
                 typename std::enable_if_t<PrelookupOption_ == lookup::prelookup::equispaced>* = nullptr>
        constexpr index_type prelookup_inside(const QueryType &q) const
        {
            //
            // Prelookup inside an equispaced grid.
            //

            const auto mm1 = _max_index - index_type{ 1u };
            return (q < _grid_vector[mm1]) ? static_cast<index_type>((q - _grid_vector[0]) / _extra_data.spacing) : mm1;

        }


        template<lookup::extrapolation ExtrapolationOption,
                 typename QueryType,
                 typename std::enable_if_t<ExtrapolationOption != lookup::extrapolation::nearest>* = nullptr>
        constexpr pure_function auto linear_frac(const QueryType &q, const prelookup_result &pr) const
        {
                const auto &idx_lo = pr[0];
                return (q - _grid_vector[idx_lo]) / (_grid_vector[idx_lo + 1] - _grid_vector[idx_lo]);

        }

        template<lookup::extrapolation ExtrapolationOption,
                 typename QueryType,
                 typename std::enable_if_t<ExtrapolationOption == lookup::extrapolation::nearest>* = nullptr>
        constexpr pure_function auto linear_frac(const QueryType &q, const prelookup_result &pr) const
        {
            const auto &e = pr[1];
            return (e == lookup::extrapolation_flag::none) ? linear_frac<lookup::extrapolation::allow>(q, pr) :
                                                             QueryType{ static_cast<typename grid_vector_type::value_type>(e) };

        }

    };


    template<class TupleOfPrelookups = std::tuple<Prelookup<>>,
             lookup::interpolation InterpolationOption = lookup::interpolation::linear,
             lookup::extrapolation ExtrapolationOption = lookup::extrapolation::allow,
             typename TableType = std::vector<double> >
    class Lookup_table
    {
        //
        // Represents a lookup table (i.e., a table of values
        // with "rank" dimensions plus a number of grid vectors)
        // in which we can interpolate.
        //

    public:

        using index_type = typename lookup::index_type;
        using table_type = TableType;
        using value_type = typename table_type::value_type;
        using tuple_of_prelookups = TupleOfPrelookups;

        static constexpr std::size_t rank = std::tuple_size_v<tuple_of_prelookups>;
        static_assert(rank != 0u, "tr7::Lookup_table: the rank must be greater than zero.");

        static constexpr lookup::interpolation interpolation_option = InterpolationOption;
        static constexpr lookup::extrapolation extrapolation_option = ExtrapolationOption;
        static_assert((interpolation_option != lookup::interpolation::nearest &&
                       interpolation_option != lookup::interpolation::previous &&
                       interpolation_option != lookup::interpolation::next) ||
                      ((interpolation_option == lookup::interpolation::nearest ||
                        interpolation_option == lookup::interpolation::previous ||
                        interpolation_option == lookup::interpolation::next) &&
                       (extrapolation_option == lookup::extrapolation::allow || extrapolation_option == lookup::extrapolation::nearest)),
                      "tr7::Lookup_table: the only extrapolation methods that are supported for the"
                      " \"nearest\", \"previous\" ant \"next\" interpolation methods are \"nearest\" or \"allow\" (and both behave the same).");


        constexpr Lookup_table() = default;

        template<typename... Args>
        Lookup_table(Args&&... args)
        {
            build(std::forward<Args>(args)...);
        }


#ifndef _MSC_VER
        template<typename Table,
                 typename... GridVectors,
                 typename std::enable_if_t<sizeof...(GridVectors) == rank>* = nullptr>
#else
        template<typename Table,
                 typename... GridVectors,
                 std::size_t Rank_ = rank,
                 typename std::enable_if_t<sizeof...(GridVectors) == Rank_>* = nullptr>
#endif
        inline void build(Table &&table_, GridVectors&&... grid_vectors)
        {
#ifdef _MSC_VER
            static_assert(Rank_ == rank, "tr7::Lookup_table::build: template parameter \"Rank_\""
                                         " should be equal to the table's rank (just allow it to get deduced!).");
#endif

            _prelookups = tuple_of_prelookups{ std::forward<GridVectors>(grid_vectors)... };
            table(std::forward<Table>(table_));

        }


        constexpr const table_type& table() const { return _table; }

        template<typename Table>
        inline void table(Table &&table)
        {
            table_impl(std::forward<Table>(table));
        }


        constexpr const tuple_of_prelookups& prelookups() const {  return _prelookups; }


        constexpr bool empty() const
        {
            return _table.empty();
        }

        constexpr void clear()
        {
            tuple_for(_prelookups, [](auto &prelookup)
                                   {
                                       if constexpr (is_tr7_Prelookup<typename std::decay_t<decltype(prelookup)> >::value) {
                                           prelookup.clear();
                                       }
                                   });
            clear_table();

        }


#ifndef _MSC_VER
        template<typename... Qs,
                 typename std::enable_if_t<sizeof...(Qs) == rank>* = nullptr>
#else
        template<typename... Qs,
                 std::size_t Rank_ = rank,
                 typename std::enable_if_t<sizeof...(Qs) == Rank_>* = nullptr>
#endif
        constexpr auto operator()(Qs&&... qs) const
        {
#ifdef _MSC_VER
            static_assert(Rank_ == rank, "tr7::Lookup_table::operator(): template parameter \"Rank_\""
                                         " should be equal to the table's rank (just allow it to get deduced!).");
#endif

            return query(std::forward<Qs>(qs)...);

        }


#ifndef _MSC_VER
        template<typename... NewGridVectors,
                 typename TableType_ = table_type,
                 typename std::enable_if_t<sizeof...(NewGridVectors) == rank &&
                                           has_resize_member_fun<TableType_>::value>* = nullptr>
#else
        template<typename... NewGridVectors,
                 typename TableType_ = table_type,
                 std::size_t Rank_ = rank,
                 typename std::enable_if_t<sizeof...(NewGridVectors) == Rank_ &&
                                           has_resize_member_fun<TableType_>::value>* = nullptr>
#endif
        inline void resample(NewGridVectors&&... new_grid_vectors)
        {
            static_assert(std::is_same_v<TableType_, table_type>, "tr7::Lookup_table::resample: template parameter \"TableType_\""
                                                                  " should be equal to the table's \"table_type\" (just allow it to get deduced!).");

#ifdef _MSC_VER
            static_assert(Rank_ == rank, "tr7::Lookup_table::resample: template parameter \"Rank_\""
                                         " should be equal to the table's rank (just allow it to get deduced!).");
#endif

            resample_table(std::forward<NewGridVectors>(new_grid_vectors)...);

        }


    private:

        tuple_of_prelookups _prelookups;
        table_type _table;
        std::array<index_type, rank - 1> _dim_strides;


        template<typename Table>
        inline void table_impl(Table &&table) noexcept(false)
        {
            //
            // When table is supplied as an array in column-major format
            // we cannot ensure a 100% safe per-dimension validation. At
            // least we'll check that the size of the table is equal to the
            // product of the sizes of the grid vectors, and hope that the
            // user has ordered the values correctly.
            //

            const auto table_size = table.size();
            std::array<std::size_t, rank> grid_vector_sizes;
            tuple_index_for(_prelookups, [&grid_vector_sizes](auto &prelookup, auto pos)
                                         {
                                             if constexpr (is_tr7_Prelookup<std::decay_t<decltype(prelookup)> >::value) {
                                                 grid_vector_sizes[pos] = prelookup.grid_vector().size();
                                             }
                                             else {
                                                 grid_vector_sizes[pos] = prelookup;
                                             }

                                         });

            const auto product_of_grid_vector_sizes = std::accumulate(std::cbegin(grid_vector_sizes), std::cend(grid_vector_sizes),
                                                                      std::size_t{ 1 }, std::multiplies<>{});

            if (product_of_grid_vector_sizes != table_size) {
                TR7_THROW_RUNTIME_ERROR("tr7::Lookup_table::validate_values: inconsistent input data. The size of"
                                        " the table is not equal to the product of the grid vectors' sizes.");

            }

            // support for the empty object
            if (product_of_grid_vector_sizes == 0u) {
                for (const auto &gvs : grid_vector_sizes) {
                    if (gvs != 0u) {
                        TR7_THROW_RUNTIME_ERROR("tr7::Lookup_table::validate_values: inconsistent input data. The table"
                                                " is empty but there's at least one non-empty grid vector.");

                    }

                }

                clear_table();
                return;

            }

            // calculate extra necessary data and assign the table
            if constexpr (rank > 1u) {
                std::inclusive_scan(std::cbegin(grid_vector_sizes), std::prev(std::cend(grid_vector_sizes)),
                                    std::begin(_dim_strides), std::multiplies<>{});

            }


            _extra_data.build(table, grid_vector_sizes, std::get<0>(_prelookups));

            if constexpr (std::is_same_v<Table, table_type>) {
                _table = std::forward<Table>(table);
            }
            else {
                _table = table_type(std::cbegin(table), std::cend(table));
            }

        }


        constexpr void clear_table()
        {
            std::fill(std::begin(_dim_strides), std::end(_dim_strides), index_type{ 0 });
            _extra_data.clear();

#ifndef _MSC_VER
            if constexpr (has_clear_member_fun<table_type>::value) {
#else
            if constexpr (is_specialization<table_type, std::vector>::value ||
                          is_specialization<table_type, std::map>::value ||
                          is_specialization<table_type, std::unordered_map>::value) {
#endif
                _table.clear();

            }
            else {
                _table = table_type{};

            }

        }


        //
        // Resampling of the table's values to a new grid impl.
        //

        template<typename NewGridVector,
                 typename TableType_ = table_type,
                 std::size_t Rank_ = rank,
                 typename std::enable_if_t<Rank_ == 1u &&
                                           has_resize_member_fun<TableType_>::value>* = nullptr>
        inline void resample_table(NewGridVector &&new_grid_vector)
        {
            //
            // Resampling for rank = 1 (super simple).
            //

            table_type new_table;
            new_table.resize(new_grid_vector.size());
            auto nt = std::begin(new_table);
            for (auto &&x : new_grid_vector) {
                *nt++ = query(x);
            }

            build(std::move(new_table), std::forward<NewGridVector>(new_grid_vector));

        }

        template<typename... NewGridVectors,
                 typename TableType_ = table_type,
                 std::size_t Rank_ = rank,
                 typename std::enable_if_t<sizeof...(NewGridVectors) == Rank_ &&
                                           Rank_ >= 2u &&
                                           has_resize_member_fun<TableType_>::value>* = nullptr>
        inline void resample_table(NewGridVectors&&... new_grid_vectors)
        {
            //
            // Resampling helper for rank >= 2: forwards the parameter pack as a tuple.
            //

            const auto ind2sub = [](const auto &dims, auto ind)
            {
                std::remove_cv_t<std::remove_reference_t<decltype(dims)> > subs{};
                auto accum{ ind };
                auto d = std::cbegin(dims);
                for (auto &&s : subs) {
                    s = accum % *d;
                    accum = (accum - s) / *d++;

                }

                return subs;

            };


            const auto new_table = [&ind2sub, this](const auto &tuple_of_new_grid_vectors)
            {
                std::array<std::size_t, rank> new_dims{};
                tuple_index_for(tuple_of_new_grid_vectors, [&new_dims](const auto &gv, auto count)
                                                           {
                                                               std::get<count>(new_dims) = gv.size();
                                                           });

                table_type new_table_;
                new_table_.resize(std::accumulate(std::cbegin(new_dims), std::cend(new_dims),
                                  std::size_t{ 1u }, std::multiplies<>{}));

                index_type count{ 0u };
                for (auto &&nt : new_table_) {
                    // calculate the sampling point
                    const auto subs = ind2sub(new_dims, count++);
                    std::tuple<typename NewGridVectors::value_type...> point{};
                    tuple_index_for(point, [&subs, &tuple_of_new_grid_vectors](auto &p, auto dim)
                                           {
                                               p = std::get<dim>(tuple_of_new_grid_vectors)[subs[dim]];
                                           });
                    nt = resampling_query(point);

                }

                return new_table_;

            };


            build(new_table(std::forward_as_tuple(new_grid_vectors...)),
                  std::forward<NewGridVectors>(new_grid_vectors)...);

        }

        template <typename TupleOfQueries, index_type... I>
        constexpr auto resampling_query(const TupleOfQueries &qs, std::integer_sequence<index_type, I...>)
        {
            return query(std::get<I>(qs)...);
        }

        template<typename TupleOfQueries>
        constexpr auto resampling_query(const TupleOfQueries &qs)
        {
            return resampling_query(qs, std::make_integer_sequence<index_type, std::tuple_size<TupleOfQueries>::value>{});
        }


        //
        // Call the prelookup impl.
        //

        template<typename PrelookupType>
        struct is_tr7_Prelookup : std::false_type {};

        template<lookup::prelookup PrelookupOption,
                 bool UsePreviousIndex, typename GridVectorType>
        struct is_tr7_Prelookup<Prelookup<PrelookupOption, UsePreviousIndex, GridVectorType>> : std::true_type {};

        template<typename QueryType, typename PrelookupType, bool RequiresLinearFrac, typename = void>
        struct is_prelookup_result : std::false_type {};

        template<typename QueryType, typename PrelookupType, bool RequiresLinearFrac>
        struct is_prelookup_result<QueryType, PrelookupType, RequiresLinearFrac,
                                   std::void_t<std::enable_if_t<(!RequiresLinearFrac && std::is_same_v<QueryType, typename PrelookupType::prelookup_result>) ||
                                                                (RequiresLinearFrac && std::is_same_v<QueryType,
                                                                                                      typename std::pair<typename PrelookupType::prelookup_result,
                                                                                                                         typename PrelookupType::grid_vector_type::value_type> >) > > > : std::true_type {};

        template<typename PrelookupType, typename QueryType,
                 typename = void>
        struct linear_frac_type {};

        template<typename PrelookupType, typename QueryType>
        struct linear_frac_type<PrelookupType, QueryType,
                                std::void_t<std::enable_if_t<!is_tr7_Prelookup<PrelookupType>::value ||
                                                             is_prelookup_result<QueryType, PrelookupType, true>::value> > >
        {
            using type = typename std::decay_t<QueryType>::second_type;
        };

        template<typename PrelookupType, typename QueryType>
        struct linear_frac_type<PrelookupType, QueryType,
                                std::void_t<std::enable_if_t<is_tr7_Prelookup<PrelookupType>::value &&
                                                             !is_prelookup_result<QueryType, PrelookupType, true>::value> > >
        {
            using type = typename std::decay_t<QueryType>;
        };

        template<typename PrelookupType, typename QueryType>
        using linear_frac_type_t = typename linear_frac_type<PrelookupType, QueryType>::type;

        template<typename... Types>
        struct tuple_of_linear_fracs {};

        template<typename ArrayOfPrelookups, typename... QueryTypes>
        struct tuple_of_linear_fracs<ArrayOfPrelookups, QueryTypes...>
        {
            using type = decay_tuple<std::tuple<linear_frac_type_t<typename ArrayOfPrelookups::value_type, QueryTypes>... > >;
        };

        template<typename... PrelookupTypes, typename... QueryTypes>
        struct tuple_of_linear_fracs<std::tuple<PrelookupTypes...>, QueryTypes...>
        {
            using type = decay_tuple<std::tuple<linear_frac_type_t<PrelookupTypes, QueryTypes>... > >;
        };

        template<bool RequiresLinearFrac,
                 typename PrelookupType,
                 typename QueryType>
        static constexpr pure_function auto evaluate_prelookup(const PrelookupType &prelookup, const QueryType &q)
        {
            static_assert((RequiresLinearFrac && (interpolation_option != lookup::interpolation::previous && interpolation_option != lookup::interpolation::next)) ||
                          (!RequiresLinearFrac && (interpolation_option == lookup::interpolation::previous || interpolation_option == lookup::interpolation::next)),
                          "tr7::Lookup_table::evaluate_prelookup: if \"RequiresLinearFrac\" is true, then"
                          " \"interpolation_option\" MUST BE \"previous\" of \"next\". Else, if \"RequiresLinearFrac\" is false"
                          " then \"interpolation_option\" CANNOT BE \"previous\" of \"next\"");

            if constexpr (is_tr7_Prelookup<PrelookupType>::value &&
                          !is_prelookup_result<QueryType, PrelookupType, RequiresLinearFrac>::value) {
                return prelookup.template evaluate<RequiresLinearFrac, extrapolation_option>(q);
            }
            else {
                return q;
            }

        }


        //
        // Linear/nearest/previous/next interpolation impl.
        //

        template<typename... QueryTypes,
                 std::size_t Rank_ = rank,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<sizeof...(QueryTypes) == Rank_ &&
                                           !is_specialization<std::common_type_t<QueryTypes...>, std::tuple>::value &&
                                           (Rank_ >= 5u ||
                                            InterpolationOption_ == lookup::interpolation::spline)>* = nullptr>
        constexpr auto query(QueryTypes&&... qs) const
        {
            //
            // Query helper for "interpN" (any N if interpolation_option = spline;
            // otherwise only used if N >= 5): forwards the parameter pack as a tuple.
            //

            return query(std::forward_as_tuple(qs...));

        }

        template<typename QueryType,
                 std::size_t Rank_ = rank,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<Rank_ == 1u &&
                                           (InterpolationOption_ == lookup::interpolation::linear ||
                                            InterpolationOption_ == lookup::interpolation::nearest ||
                                            InterpolationOption_ == lookup::interpolation::previous ||
                                            InterpolationOption_ == lookup::interpolation::next)>* = nullptr>
        constexpr auto query(const QueryType &xq) const
        {
            //
            // "interp1" with interpolation_option = linear/nearest/previous/next.
            //

            if constexpr (InterpolationOption_ == lookup::interpolation::linear ||
                          InterpolationOption_ == lookup::interpolation::nearest) {

                const auto [pr, frac_x] = evaluate_prelookup<true>(std::get<0>(_prelookups), xq);
                const auto &i_lo = pr[0];

                if constexpr (InterpolationOption_ == lookup::interpolation::linear) {
                    return line_formula(_table[i_lo + 1], _table[i_lo], frac_x);

                }
                else {
                    return frac_x > decltype(frac_x){ 0.5 } ? _table[i_lo + 1u] : _table[i_lo];

                }

            }
            else {
                const auto pr = evaluate_prelookup<false>(std::get<0>(_prelookups), xq);
                const auto &i_lo = pr[0];
                const auto &e_x = pr[1];

                if constexpr (InterpolationOption_ == lookup::interpolation::previous) {
                    return e_x == lookup::extrapolation_flag::upper ? _table[i_lo + 1u] : _table[i_lo];

                }
                else {
                    return e_x != lookup::extrapolation_flag::lower ? _table[i_lo + 1u] : _table[i_lo];

                }

            }

        }

        template<typename XQueryType, typename YQueryType,
                 std::size_t Rank_ = rank,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<Rank_ == 2u &&
                                           (InterpolationOption_ == lookup::interpolation::linear ||
                                            InterpolationOption_ == lookup::interpolation::nearest ||
                                            InterpolationOption_ == lookup::interpolation::previous ||
                                            InterpolationOption_ == lookup::interpolation::next)>* = nullptr>
        constexpr auto query(const XQueryType &xq, const YQueryType &yq) const
        {
            //
            // "interp2" with interpolation_option = linear/nearest/previous/next.
            //

            if constexpr (InterpolationOption_ == lookup::interpolation::linear ||
                          InterpolationOption_ == lookup::interpolation::nearest) {
                auto [pr_x, frac_x] = evaluate_prelookup<true>(std::get<0>(_prelookups), xq);
                auto &i_lo = pr_x[0];
                auto [pr_y, frac_y] = evaluate_prelookup<true>(std::get<1>(_prelookups), yq);
                auto &j_lo = pr_y[0];

                if constexpr (InterpolationOption_ == lookup::interpolation::linear) {
                    auto ij_lo = i_lo + _dim_strides[0] * j_lo;

                    const auto fL = line_formula(_table[ij_lo + 1], _table[ij_lo], frac_x);

                    ij_lo += _dim_strides[0];

                    const auto fR = line_formula(_table[ij_lo + 1], _table[ij_lo], frac_x);

                    return line_formula(fR, fL, frac_y);

                }
                else {
                    if (frac_x > decltype(frac_x){ 0.5 }) {
                        ++i_lo;
                    }
                    if (frac_y > decltype(frac_y){ 0.5 }) {
                        ++j_lo;
                    }

                    return _table[i_lo + _dim_strides[0] * j_lo];

                }

            }
            else {
                auto pr_x = evaluate_prelookup<false>(std::get<0>(_prelookups), xq);
                auto &i_lo = pr_x[0];
                const auto &e_x = pr_x[1];
                auto pr_y = evaluate_prelookup<false>(std::get<1>(_prelookups), yq);
                auto &j_lo = pr_y[0];
                const auto &e_y = pr_y[1];

                if constexpr (InterpolationOption_ == lookup::interpolation::previous) {
                    if (e_x == lookup::extrapolation_flag::upper) {
                        ++i_lo;
                    }
                    if (e_y == lookup::extrapolation_flag::upper) {
                        ++j_lo;
                    }

                }
                else {
                    if (e_x != lookup::extrapolation_flag::lower) {
                        ++i_lo;
                    }
                    if (e_y != lookup::extrapolation_flag::lower) {
                        ++j_lo;
                    }

                }

                return _table[i_lo + _dim_strides[0] * j_lo];

            }

       }

        template<typename XQueryType, typename YQueryType, typename ZQueryType,
                 std::size_t Rank_ = rank,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<Rank_ == 3u &&
                                           (InterpolationOption_ == lookup::interpolation::linear ||
                                            InterpolationOption_ == lookup::interpolation::nearest ||
                                            InterpolationOption_ == lookup::interpolation::previous ||
                                            InterpolationOption_ == lookup::interpolation::next)>* = nullptr>
        constexpr auto query(const XQueryType &xq, const YQueryType &yq, const ZQueryType &zq) const
        {
            //
            // "interp3" with interpolation_option = linear/nearest/previous/next.
            //

            if constexpr (InterpolationOption_ == lookup::interpolation::linear ||
                          InterpolationOption_ == lookup::interpolation::nearest) {

                auto [pr_x, frac_x] = evaluate_prelookup<true>(std::get<0>(_prelookups), xq);
                auto &i_lo = pr_x[0];
                auto [pr_y, frac_y] = evaluate_prelookup<true>(std::get<1>(_prelookups), yq);
                auto &j_lo = pr_y[0];
                auto [pr_z, frac_z] = evaluate_prelookup<true>(std::get<2>(_prelookups), zq);
                auto &k_lo = pr_z[0];

                if constexpr (InterpolationOption_ == lookup::interpolation::linear) {
                    auto ijk_lo = i_lo + _dim_strides[0] * j_lo + _dim_strides[1] * k_lo;

                    const auto fL = line_formula(line_formula(_table[ijk_lo + 1 + _dim_strides[0]], _table[ijk_lo + _dim_strides[0]], frac_x),
                                                 line_formula(_table[ijk_lo + 1], _table[ijk_lo], frac_x),
                                                 frac_y);

                    ijk_lo += _dim_strides[1];

                    const auto fR = line_formula(line_formula(_table[ijk_lo + 1 + _dim_strides[0]], _table[ijk_lo + _dim_strides[0]], frac_x),
                                                 line_formula(_table[ijk_lo + 1], _table[ijk_lo], frac_x),
                                                 frac_y);

                    return line_formula(fR, fL, frac_z);

                }
                else {
                    if (frac_x > decltype(frac_x){ 0.5 }) {
                        ++i_lo;
                    }
                    if (frac_y > decltype(frac_y){ 0.5 }) {
                        ++j_lo;
                    }
                    if (frac_z > decltype(frac_z){ 0.5 }) {
                        ++k_lo;
                    }

                    return _table[i_lo + _dim_strides[0] * j_lo + _dim_strides[1] * k_lo];

                }

            }
            else {
                auto pr_x = evaluate_prelookup<false>(std::get<0>(_prelookups), xq);
                auto &i_lo = pr_x[0];
                const auto &e_x = pr_x[1];
                auto pr_y = evaluate_prelookup<false>(std::get<1>(_prelookups), yq);
                auto &j_lo = pr_y[0];
                const auto &e_y = pr_y[1];
                auto pr_z = evaluate_prelookup<false>(std::get<2>(_prelookups), zq);
                auto &k_lo = pr_z[0];
                const auto &e_z = pr_z[1];

                if constexpr (InterpolationOption_ == lookup::interpolation::previous) {
                    if (e_x == lookup::extrapolation_flag::upper) {
                        ++i_lo;
                    }
                    if (e_y == lookup::extrapolation_flag::upper) {
                        ++j_lo;
                    }
                    if (e_z == lookup::extrapolation_flag::upper) {
                        ++k_lo;
                    }

                }
                else {
                    if (e_x != lookup::extrapolation_flag::lower) {
                        ++i_lo;
                    }
                    if (e_y != lookup::extrapolation_flag::lower) {
                        ++j_lo;
                    }
                    if (e_z != lookup::extrapolation_flag::lower) {
                        ++k_lo;
                    }

                }

                return _table[i_lo + _dim_strides[0] * j_lo + _dim_strides[1] * k_lo];

            }

        }

        template<typename XQueryType, typename YQueryType, typename ZQueryType, typename TQueryType,
                 std::size_t Rank_ = rank,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<Rank_ == 4u &&
                                           (InterpolationOption_ == lookup::interpolation::linear ||
                                            InterpolationOption_ == lookup::interpolation::nearest ||
                                            InterpolationOption_ == lookup::interpolation::previous ||
                                            InterpolationOption_ == lookup::interpolation::next)>* = nullptr>
        constexpr auto query(const XQueryType &xq, const YQueryType &yq, const ZQueryType &zq, const TQueryType &tq) const
        {
            //
            // "interp4" with interpolation_option = linear/nearest/previous/next.
            //

            if constexpr (InterpolationOption_ == lookup::interpolation::linear ||
                          InterpolationOption_ == lookup::interpolation::nearest) {

                auto [pr_x, frac_x] = evaluate_prelookup<true>(std::get<0>(_prelookups), xq);
                auto &i_lo = pr_x[0];
                auto [pr_y, frac_y] = evaluate_prelookup<true>(std::get<1>(_prelookups), yq);
                auto &j_lo = pr_y[0];
                auto [pr_z, frac_z] = evaluate_prelookup<true>(std::get<2>(_prelookups), zq);
                auto &k_lo = pr_z[0];
                auto [pr_t, frac_t] = evaluate_prelookup<true>(std::get<3>(_prelookups), tq);
                auto &t_lo = pr_t[0];

                if constexpr (InterpolationOption_ == lookup::interpolation::linear) {
                    auto ijkt_lo = i_lo + _dim_strides[0] * j_lo + _dim_strides[1] * k_lo + _dim_strides[2] * t_lo;

                    auto mL = line_formula(line_formula(_table[ijkt_lo + 1 + _dim_strides[0]], _table[ijkt_lo + _dim_strides[0]], frac_x),
                                           line_formula(_table[ijkt_lo + 1], _table[ijkt_lo], frac_x),
                                           frac_y);

                    ijkt_lo += _dim_strides[1];

                    auto fR = line_formula(line_formula(_table[ijkt_lo + 1 + _dim_strides[0]], _table[ijkt_lo + _dim_strides[0]], frac_x),
                                           line_formula(_table[ijkt_lo + 1], _table[ijkt_lo], frac_x),
                                           frac_y);

                    mL = line_formula(fR, mL, frac_z);

                    ijkt_lo += _dim_strides[2];

                    auto mR = line_formula(line_formula(_table[ijkt_lo + 1 + _dim_strides[0]], _table[ijkt_lo + _dim_strides[0]], frac_x),
                                           line_formula(_table[ijkt_lo + 1], _table[ijkt_lo], frac_x),
                                           frac_y);

                    ijkt_lo -= _dim_strides[1];

                    fR = line_formula(line_formula(_table[ijkt_lo + 1 + _dim_strides[0]], _table[ijkt_lo + _dim_strides[0]], frac_x),
                                      line_formula(_table[ijkt_lo + 1], _table[ijkt_lo], frac_x),
                                      frac_y);

                    mR = line_formula(mR, fR, frac_z);

                    return line_formula(mR, mL, frac_t);

                }
                else {
                    if (frac_x > decltype(frac_x){ 0.5 }) {
                        ++i_lo;
                    }
                    if (frac_y > decltype(frac_y){ 0.5 }) {
                        ++j_lo;
                    }
                    if (frac_z > decltype(frac_z){ 0.5 }) {
                        ++k_lo;
                    }
                    if (frac_t > decltype(frac_t){ 0.5 }) {
                        ++t_lo;
                    }

                    return _table[i_lo + _dim_strides[0] * j_lo + _dim_strides[1] * k_lo + _dim_strides[2] * t_lo];

                }

            }
            else {
                auto pr_x = evaluate_prelookup<false>(std::get<0>(_prelookups), xq);
                auto &i_lo = pr_x[0];
                const auto &e_x = pr_x[1];
                auto pr_y = evaluate_prelookup<false>(std::get<1>(_prelookups), yq);
                auto &j_lo = pr_y[0];
                const auto &e_y = pr_y[1];
                auto pr_z = evaluate_prelookup<false>(std::get<2>(_prelookups), zq);
                auto &k_lo = pr_z[0];
                const auto &e_z = pr_z[1];
                auto pr_t = evaluate_prelookup<false>(std::get<3>(_prelookups), tq);
                auto &t_lo = pr_t[0];
                const auto &e_t = pr_t[1];

                if constexpr (InterpolationOption_ == lookup::interpolation::previous) {
                    if (e_x == lookup::extrapolation_flag::upper) {
                        ++i_lo;
                    }
                    if (e_y == lookup::extrapolation_flag::upper) {
                        ++j_lo;
                    }
                    if (e_z == lookup::extrapolation_flag::upper) {
                        ++k_lo;
                    }
                    if (e_t == lookup::extrapolation_flag::upper) {
                        ++t_lo;
                    }

                }
                else {
                    if (e_x != lookup::extrapolation_flag::lower) {
                        ++i_lo;
                    }
                    if (e_y != lookup::extrapolation_flag::lower) {
                        ++j_lo;
                    }
                    if (e_z != lookup::extrapolation_flag::lower) {
                        ++k_lo;
                    }
                    if (e_t != lookup::extrapolation_flag::lower) {
                        ++t_lo;
                    }

                }

                return _table[i_lo + _dim_strides[0] * j_lo + _dim_strides[1] * k_lo + _dim_strides[2] * t_lo];

            }

        }

        template<typename... QueryTypes,
                 std::size_t Rank_ = rank,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<Rank_ >= 5u &&
                                           sizeof... (QueryTypes) == Rank_ &&
                                           (InterpolationOption_ == lookup::interpolation::linear ||
                                            InterpolationOption_ == lookup::interpolation::nearest ||
                                            InterpolationOption_ == lookup::interpolation::previous ||
                                            InterpolationOption_ == lookup::interpolation::next)>* = nullptr>
        constexpr auto query(const std::tuple<QueryTypes...> &qs) const
        {
            //
            // "interpN" (N >= 5) with interpolation_option = linear/nearest/previous/next.
            //

            if constexpr (InterpolationOption_ == lookup::interpolation::linear) {
                const auto idx_lo_and_frac_at_pos = [this](auto &frac, const auto &prelookup, const auto &q)
                                                    {
                                                        const auto [pr, frac_] = evaluate_prelookup<true>(prelookup, q);
                                                        frac = frac_;
                                                        return pr[0];
                                                    
                                                    };

                typename tuple_of_linear_fracs<tuple_of_prelookups, QueryTypes...>::type fracs{};
                auto idx = idx_lo_and_frac_at_pos(std::get<0>(fracs), std::get<0>(_prelookups), std::get<0>(qs));
                index_for<rank - 1u>([&idx, &qs, &fracs, idx_lo_and_frac_at_pos, this](auto pos)
                                     {
                                         idx += idx_lo_and_frac_at_pos(std::get<pos + 1u>(fracs),
                                                                       std::get<pos + 1u>(_prelookups),
                                                                       std::get<pos + 1u>(qs)) * _dim_strides[pos];

                                     });

#ifdef TR7_WITH_AD
                // with AD, the lookup table's output should be an AD-type, not a
                // "Lookup_table::value_type" (i.e., the data that the table holds,
                // which isn't an AD-type). For now, we'll only support tables
                // holding scalars or std::array values...
                if constexpr (tr7::any_tuple_type_is_AD_scalar<decltype(fracs)>::value) {
                      if constexpr (is_std_array<value_type>::value) {
                            std::array<typename tr7::tuple_common_type<decltype(fracs)>::type,
                                       std::tuple_size_v<value_type> > ret{};
                            recursive_linear_queryN<rank>(ret, idx, fracs);
                            return ret;

                        }
                        else {
                            typename tr7::tuple_common_type<decltype(fracs)>::type ret{};
                            recursive_linear_queryN<rank>(ret, idx, fracs);
                            return ret;

                        }

                }
                else {
#endif
                    value_type ret{};
                    recursive_linear_queryN<rank>(ret, idx, fracs);
                    return ret;

#ifdef TR7_WITH_AD
                }
#endif

            }
            else {
                const auto idx_lo_at_pos = [this](const auto &prelookup, const auto &q)
                {
                    //auto [idx_lo, e] = prelookup(q);
                    if constexpr (InterpolationOption_ == lookup::interpolation::nearest) {
                        auto [pr, frac] = evaluate_prelookup<true>(prelookup, q);
                        auto &idx_lo = pr[0];

                        if (frac > decltype(frac){ 0.5 }) {
                            ++idx_lo;
                        }

                        return idx_lo;

                    }
                    else {
                        auto pr = evaluate_prelookup<false>(prelookup, q);
                        auto &idx_lo = pr[0];
                        const auto &e = pr[1];

                        if constexpr (InterpolationOption_ == lookup::interpolation::previous) {
                            if (e == lookup::extrapolation_flag::upper) {
                                ++idx_lo;
                            }

                        }
                        else {
                            if (e != lookup::extrapolation_flag::lower) {
                                ++idx_lo;
                            }

                        }

                        return idx_lo;

                    }

                };

                auto idx = idx_lo_at_pos(std::get<0>(_prelookups), std::get<0>(qs));
                index_for<rank - 1u>([&idx, &qs, idx_lo_at_pos, this](auto pos)
                                     {
                                         idx += idx_lo_at_pos(std::get<pos + 1u>(_prelookups),
                                                              std::get<pos + 1u>(qs)) * _dim_strides[pos];

                                     });

                return _table[idx];

            }

        }

        template<std::size_t Count,
                 typename ValueType,
                 typename TupleOfFracs,
                 std::size_t Rank_ = rank,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<Rank_ >= 5u &&
                                           Count >= 1u &&
                                           std::tuple_size_v<TupleOfFracs> == Rank_ &&
                                           InterpolationOption_ == lookup::interpolation::linear>* = nullptr>
        constexpr void recursive_linear_queryN(ValueType &f, index_type idx_lo, const TupleOfFracs &fracs) const
        {
            //
            // A helper function for "interpN" (N >= 5) with
            // interpolation_option = linear.
            //

            ValueType fR{};
            if constexpr (Count - 1u != 0u) {
                recursive_linear_queryN<Count - 1u>(f, idx_lo, fracs);
                recursive_linear_queryN<Count - 1u>(fR, idx_lo + std::get<Count - 2u>(_dim_strides), fracs);

            }
            else {
                f = _table[idx_lo];
                fR = _table[idx_lo + 1];

            }

            f = line_formula(fR, f, std::get<Count - 1u>(fracs));

        };


        //
        // Linear interpolation support for iterable tables that don't define
        // the adequate operators (I'm thinking about suitable types to represent
        // a "value" in a table, e.g. std::array, std:vector, std::map, and
        // std::unordered_map).
        //

        template<typename ValueType, typename = void>
        struct has_minus_op : std::false_type {};

        template<typename ValueType>
        struct has_minus_op<ValueType, std::void_t<decltype(std::declval<ValueType>() - std::declval<ValueType>())> > : std::true_type {};

        template<typename ValueType, typename = void>
        struct has_plus_op : std::false_type {};

        template<typename ValueType>
        struct has_plus_op<ValueType, std::void_t<decltype(std::declval<ValueType>() + std::declval<ValueType>())> > : std::true_type {};

        template<typename ValueType, typename Scalar, typename = void>
        struct has_timesscalar_op : std::false_type {};

        template<typename ValueType, typename Scalar>
        struct has_timesscalar_op<ValueType, Scalar, std::void_t<decltype(std::declval<Scalar>() * std::declval<ValueType>())> > : std::true_type {};

        template<typename ValueType, typename Scalar, typename = void>
        struct has_linear_interp_ops : std::false_type {};

        template<typename ValueType, typename Scalar>
        struct has_linear_interp_ops<ValueType, Scalar, std::void_t<std::enable_if_t<has_minus_op<ValueType>::value &&
                                                                                     has_plus_op<ValueType>::value &&
                                                                                     has_timesscalar_op<ValueType, Scalar>::value> > > : std::true_type {};

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_linear_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::linear>* = nullptr>
        static constexpr pure_function auto line_formula(const ValueType &b, const ValueType &a, const Scalar &f)
        {
            return a + f * (b - a);
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_linear_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::linear &&
                                           !(is_specialization<ValueType, std::map>::value ||
                                             is_specialization<ValueType, std::unordered_map>::value)>* = nullptr>
        static constexpr pure_function auto line_formula(const ValueType &b, ValueType a, const Scalar &f)
        {
            auto bi = std::cbegin(b);
            for (auto &&ai : a) {
                ai += f * (*bi++ - ai);
            }
            return a;

        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_linear_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::linear &&
                                           (is_specialization<ValueType, std::map>::value ||
                                            is_specialization<ValueType, std::unordered_map>::value)>* = nullptr>
        static constexpr pure_function auto line_formula(const ValueType &b, ValueType a, const Scalar &f)
        {
            // std::map and std::unordered_map have std::pair iterators...
            auto bi = std::cbegin(b);
            for (auto &&ai : a) {
                ai.second += f * (bi->second - ai.second);
                ++bi;

            }

            return a;

        }


#ifdef TR7_WITH_AD
        // with AD, the lookup table's output should be of "ScalarType" (which is an AD-scalar,
        // since it comes from a query), not "ValueType" (i.e., the data that the table holds,
        // which isn't an AD-type). For now, we'll only support tables holding std::array values...
        template<typename ValueType,
                 typename Scalar,
                 std::size_t N,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_linear_interp_ops<ValueType, Scalar>::value &&
                                           tr7::is_AD_scalar_v<typename std::decay_t<Scalar> > &&
                                           InterpolationOption_ == lookup::interpolation::linear>* = nullptr>
        static constexpr pure_function auto line_formula(const std::array<ValueType, N> &b, const std::array<ValueType, N> &a, const Scalar &f)
        {
            std::array<Scalar, N> ret{};
            tr7::tuple_index_for(ret, [&a, &b, &f](auto &ret_i, auto i)
                                      {
                                          ret_i = a[i] + f * (b[i] - a[i]);
                                      });
            return ret;

        }

#endif


        //
        // Cubic spline interpolation impl., for all ranks. This method is already expensive enough,
        // so we'll not bother creating taylored versions for rank <= 4, unlike we did with the
        // other ones. It will be an "interpN" for all N's (Simulink does it like this).
        //

        template<typename... QueryTypes,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        constexpr auto query(const std::tuple<QueryTypes...> &qs) const
        {
            //
            // Natural cubic spline interpolation, plagiarized from Simulink's "Lookup Table" block.
            // NOTE: Matlab's built-in "griddedInterpolant" is MUCH FASTER than this (x 10), I suppose
            // it is because it caches the interpolation weights it calculates or something like that...
            //

            bool buf_bank = false;
            auto *yy = &_extra_data.yyA_work[0];
            auto *yy2 = &_extra_data.yy2_work[0];
            const auto *v = &_table[0];
            const auto *y2 = &_extra_data.y20_container[0];
            tuple_index_for(_prelookups, [&qs, &buf_bank, &v, &y2, &yy, &yy2, this](const auto &prelookup, auto dim)
                                         {
                                             const auto &q = std::get<dim>(qs);
                                             const auto &gv = prelookup.grid_vector();
                                             const auto [pr, p] = evaluate_prelookup<true>(prelookup, q);
                                             const auto &il = pr[0];

                                             using grid_type = typename std::decay_t<decltype(gv)>::value_type;
                                             const auto len = prelookup.max_index() + 1u;
                                             const auto iu = il + 1u;
                                             const auto h = gv[iu] - gv[il];
                                             const auto h2_div6 = h * h / grid_type{ 6 };
                                             const auto s = grid_type{ 1 } - p;
                                             const auto pmsq = p * (p * p - grid_type{ 1 });
                                             const auto smsq = s * (s * s - grid_type{ 1 });

                                             // calculate spline curves for input in this
                                             // dimension at each value of the higher
                                             // other dimensions' points in the table
                                             if constexpr (extrapolation_option == lookup::extrapolation::linear) {
                                                 if (p > grid_type{ 1 }) {
                                                     for (auto i = 0u; i < std::get<dim>(_extra_data.elts); ++i) {
                                                         // = v[iu] - v[il] + (y2[il] * h * h) * (1.0 / 6.0); (cannot be const auto!)
                                                         auto slope = slope_formula(v[iu], v[il], h2_div6, y2[il]);

                                                         // = v[iu] + (p - 1.0) * slope;
                                                         yy[i] = yy_linextrap_formula(v[iu], p - grid_type{ 1 }, slope);

                                                         v += len;
                                                         y2 += len;
                                                     }

                                                 }
                                                 else if (p < grid_type{ 0 }) {
                                                     for (auto i = 0u; i < std::get<dim>(_extra_data.elts); ++i) {
                                                         // = v[iu] - v[il] - (y2[iu] * h * h) * (1.0 / 6.0); (cannot be const auto!)
                                                         auto slope = slope_formula(v[iu], v[il], -h2_div6, y2[iu]);

                                                         // = v[il] + p * slope;
                                                         yy[i] = yy_linextrap_formula(v[il], p, slope);

                                                         v += len;
                                                         y2 += len;
                                                     }

                                                 }
                                                 else {
                                                     for (auto i = 0u; i < std::get<dim>(_extra_data.elts); ++i) {
                                                         // = v[il] + p * (v[iu] - v[il]) + ((smsq * y2[il] + pmsq * y2[iu]) * h * h) * (1.0 / 6.0);
                                                         yy[i] = yy_formula(v[il], p, v[iu], smsq, y2[il], pmsq, y2[iu], h2_div6);

                                                         v += len;
                                                         y2 += len;
                                                     }

                                                 }

                                             }
                                             else {
                                                 for (auto i = 0u; i < std::get<dim>(_extra_data.elts); ++i) {
                                                     // = v[il] + p * (v[iu] - v[il]) + ((smsq * y2[il] + pmsq * y2[iu]) * h * h) * (1.0 / 6.0);
                                                     yy[i] = yy_formula(v[il], p, v[iu], smsq, y2[il], pmsq, y2[iu], h2_div6);

                                                     v += len;
                                                     y2 += len;
                                                 }

                                             }

                                             // set pointers to new result and calculate second derivatives
                                             if constexpr (dim < rank - 1u) {
                                                 auto *v_ = yy;
                                                 auto *y2_ = yy2;
                                                 const auto &prelookup_next = std::get<dim + 1u>(_prelookups);
                                                 const auto &max_index_next = prelookup_next.max_index();
                                                 const auto len_next = max_index_next + 1u;
                                                 const auto &gv_next = prelookup_next.grid_vector();
                                                 for (auto i = 0u; i < std::get<dim + 1u>(_extra_data.elts); ++i) {
                                                     cubic_spline_second_derivatives(gv_next, max_index_next, v_, y2_, _extra_data.u_work, _extra_data.zero);
                                                     v_ += len_next;
                                                     y2_ += len_next;

                                                 }

                                             }

                                             // set work vectors v, y2 and yy for next iteration;
                                             // the yy just calculated becomes the v in the
                                             // next iteration, y2 was just calculated for these
                                             // new points and the yy buffer is swapped to the space
                                             // for storing the next iteration's results
                                             v = yy;
                                             y2 = yy2;

                                             // swap buffers for next dimension and
                                             // toggle buf_bank for next iteration
                                             yy = buf_bank ? &_extra_data.yyB_work[0] : &_extra_data.yyA_work[0];
                                             buf_bank = !buf_bank;

                                         });

            return *v;

        }


        template<typename GridVector,
                 typename Itv, typename It2,
                 typename WorkVector, typename Z,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function void cubic_spline_second_derivatives(const GridVector &gv, index_type max_index,
                                                                            Itv vals, It2 y2, WorkVector &u, const Z &zero_)
        {
            // back boundary condition: natural spline (y''[0] == 0)
            *y2 = zero_;
            u[0] = zero_;

            // an alternative to the natural back boundary condition is specifying
            // a prescribed value for the first derivative:
            // *y2 = value_type{ -0.5 };
            // u[0] = (3.0 / (gv[1] - gv[0])) * ((vals[1] - vals[0]) / (gv[1] - gv[0]) - PRESCRIBED_FIRST_DERIVATIVE_VALUE_BACK);

            // forward substitution (tridiagonal algorithm)
            for (auto i = 1u; i < max_index; ++i) {
                const auto dxm1 = gv[i] - gv[i - 1];
                const auto dxp1 = gv[i + 1] - gv[i];
                const auto dxpm = dxp1 + dxm1;
                const auto sig = dxm1 / dxpm;

                // = sig * y2[i - 1] + 2.;
                const auto p = p_formula(sig, y2[i - 1]);

                // = (sig - 1.0) / p;
                y2[i] = y2_forward_formula(sig, p);

                // = (vals[i + 1] - vals[i]) / dxp1 - (vals[i] - vals[i - 1]) / dxm1; (cannot be const auto!)
                auto utemp = utemp_formula(vals[i + 1], vals[i], dxp1, vals[i - 1], dxm1);

                // = (((6.0 * utemp) / dxpm) - (sig * u[i - 1])) / p;
                u[i] = u_forward_formula(utemp, dxpm, sig, u[i - 1], p);

            }

            // front boundary condition: natural spline (y''[end] == 0)
            y2[max_index] = zero_;

            // an alternative to the natural front boundary condition is specifying
            // a prescribed value for the first derivative:
            // static constexpr auto qn = value_type{ 0.5 };
            // const auto un = (3.0 / (gv[max_index] - gv[max_index - 1])) * (PRESCRIBED_FIRST_DERIVATIVE_VALUE_FRONT - (vals[max_index] - vals[max_index - 1]) / (gv[max_index] - gv[max_index - 1]));
            // y2[max_index] = (un - qn * u[max_index - 1]) / (qn * y2[max_index - 1] + 1.0);

            // backward substitution (tridiagonal algorithm)
            for (auto i = max_index; i > 0u; --i) {
                // set y2[i - 1] = y2[i - 1] * y2[i] + u[i - 1];
                y2_backward_subroutine(y2[i - 1], y2[i], u[i - 1]);

            }

        }


        //
        // Cubic spline interpolation support for iterable value types that don't
        // define the adequate operators (I'm thinking about suitable types to represent
        // a "value" in a table, e.g. std::array, std:vector, std::map, and
        // std::unordered_map).
        // 

        template<typename ValueType, typename Scalar, typename = void>
        struct has_divscalar_op : std::false_type {};

        template<typename ValueType, typename Scalar>
        struct has_divscalar_op<ValueType, Scalar, std::void_t<decltype(std::declval<ValueType>() / std::declval<Scalar>())> > : std::true_type {};

        template<typename ValueType, typename Scalar, typename = void>
        struct has_scalardiv_op : std::false_type {};

        template<typename ValueType, typename Scalar>
        struct has_scalardiv_op<ValueType, Scalar, std::void_t<decltype(std::declval<Scalar>() / std::declval<ValueType>())> > : std::true_type {};

        template<typename ValueType, typename = void>
        struct has_unary_plus : std::false_type {};

        template<typename ValueType>
        struct has_unary_plus<ValueType, std::void_t<decltype(std::declval<ValueType &>() += std::declval<ValueType>())> > : std::true_type {};

        template<typename ValueType, typename = void>
        struct has_unary_prod : std::false_type {};

        template<typename ValueType>
        struct has_unary_prod<ValueType, std::void_t<decltype(std::declval<ValueType &>() *= std::declval<ValueType>())> > : std::true_type {};

        template<typename ValueType, typename = void>
        struct has_div_op : std::false_type {};

        template<typename ValueType>
        struct has_div_op<ValueType, std::void_t<decltype(std::declval<ValueType &&>() / std::declval<ValueType>())> > : std::true_type {};

        template<typename ValueType, typename Scalar, typename = void>
        struct has_unary_timesscalar_prod : std::false_type {};

        template<typename ValueType, typename Scalar>
        struct has_unary_timesscalar_prod<ValueType, Scalar, std::void_t<decltype(std::declval<ValueType &>() *= std::declval<Scalar>())> > : std::true_type {};

        template<typename ValueType, typename Scalar, typename = void>
        struct has_spline_interp_ops : std::false_type {};

        template<typename ValueType, typename Scalar>
        struct has_spline_interp_ops<ValueType, Scalar, std::void_t<std::enable_if_t<has_linear_interp_ops<ValueType, Scalar>::value &&
                                                                                     has_divscalar_op<ValueType, Scalar>::value &&
                                                                                     has_scalardiv_op<ValueType, Scalar>::value &&
                                                                                     has_unary_plus<ValueType>::value &&
                                                                                     has_unary_prod<ValueType>::value &&
                                                                                     has_div_op<ValueType>::value &&
                                                                                     has_unary_timesscalar_prod<ValueType, Scalar>::value> > > : std::true_type {};


        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_spline_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto p_formula(const Scalar &sig, const ValueType &y2_im1)
        {
            // p = sig * y2[i - 1] + 2.;
            return sig * y2_im1 + ValueType{ 2 };
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto p_formula(const Scalar &sig, ValueType y2_im1)
        {
            // p = sig * y2[i - 1] + 2.;
            for (auto &&vi : y2_im1) {
                (vi *= sig) += decltype(vi){ 2 };
            }
            return y2_im1;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto p_formula(const Scalar &sig, ValueType y2_im1)
        {
            // p = sig * y2[i - 1] + 2.;
            for (auto &&vi : y2_im1) {
                (vi.second *= sig) += decltype(vi.second){ 2 };
            }
            return y2_im1;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_spline_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto y2_forward_formula(const Scalar &sig, const ValueType &p)
        {
            // y2 = (sig - 1.) / p;
            return (sig - Scalar{ 1 }) / p;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto y2_forward_formula(const Scalar &sig, ValueType p)
        {
            // y2 = (sig - 1.) / p;
            const auto sm1 = sig - Scalar{ 1 };
            for (auto &&vi : p) {
                vi = sm1 / vi;
            }
            return p;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto y2_forward_formula(const Scalar &sig, ValueType p)
        {
            // y2 = (sig - 1.) / p;
            const auto sm1 = sig - Scalar{ 1 };
            for (auto &&vi : p) {
                vi.second = sm1 / vi.second;
            }
            return p;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_spline_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto utemp_formula(const ValueType &vals_ip1, const ValueType &vals_i,
                                                               const Scalar &dxp1, const ValueType &vals_im1, const Scalar &dxm1)
        {
            // utemp = (vals[i + 1] - vals[i]) / dxp1 - (vals[i] - vals[i - 1]) / dxm1;
            return (vals_ip1 - vals_i) / dxp1 - (vals_i - vals_im1) / dxm1;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto utemp_formula(const ValueType &vals_ip1, ValueType vals_i,
                                                          const Scalar &dxp1, const ValueType &vals_im1, const Scalar &dxm1)
        {
            // utemp = (vals[i + 1] - vals[i]) / dxp1 - (vals[i] - vals[i - 1]) / dxm1;
            auto v_ip1 = std::cbegin(vals_ip1);
            auto v_im1 = std::cbegin(vals_im1);
            for (auto &&v_i : vals_i) {
                v_i = (*v_ip1++ - v_i) / dxp1 - (v_i - *v_im1++) / dxm1;
            }
            return vals_i;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto utemp_formula(const ValueType &vals_ip1, ValueType vals_i,
                                                          const Scalar &dxp1, const ValueType &vals_im1, const Scalar &dxm1)
        {
            // utemp = (vals[i + 1] - vals[i]) / dxp1 - (vals[i] - vals[i - 1]) / dxm1;
            auto v_ip1 = std::cbegin(vals_ip1);
            auto v_im1 = std::cbegin(vals_im1);
            for (auto &&v_i : vals_i) {
                v_i.second = (v_ip1->second - v_i.second) / dxp1 - (v_i.second - v_im1->second) / dxm1;
                ++v_ip1;
                ++v_im1;
            }
            return vals_i;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_spline_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto u_forward_formula(const ValueType &utemp, const Scalar &dxpm,
                                                              const Scalar &sig, const ValueType &u_im1, const ValueType &p)
        {
            // u[i] = (((6.0 * utemp) / dxpm) - (sig * u[i - 1])) / p;
            return (Scalar{ 6 } * utemp / dxpm - (sig * u_im1)) / p;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto& u_forward_formula(ValueType &utemp, const Scalar &dxpm,
                                                               const Scalar &sig, const ValueType &u_im1_, const ValueType &p_)
        {
            // u[i] = (((6.0 * utemp) / dxpm) - (sig * u[i - 1])) / p;
            const auto six = Scalar{ 6 } / dxpm;
            auto u_im1 = std::cbegin(u_im1_);
            auto p = std::cbegin(p_);
            for (auto &&u : utemp) {
                u = (six * u - sig * (*u_im1++)) / *p++;
            }
            return utemp;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto& u_forward_formula(ValueType &utemp, const Scalar &dxpm,
                                                               const Scalar &sig, const ValueType &u_im1_, const ValueType &p_)
        {
            // u[i] = (((6.0 * utemp) / dxpm) - (sig * u[i - 1])) / p;
            const auto six = Scalar{ 6 } / dxpm;
            auto u_im1 = std::cbegin(u_im1_);
            auto p = std::cbegin(p_);
            for (auto &&u : utemp) {
                u.second = (six * u.second - sig * (u_im1->second)) / p->second;
                ++p; 
                ++u_im1;
            }
            return utemp;
        }

        template<typename ValueType,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<(has_unary_plus<ValueType>::value && has_unary_prod<ValueType>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function void y2_backward_subroutine(ValueType &y2_im1, const ValueType &y2_i, const ValueType &u_im1)
        {
            // y2[i - 1] = y2[i - 1] * y2[i] + u[i - 1];
            (y2_im1 *= y2_i) += u_im1;
        }

        template<typename ValueType,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!(has_unary_plus<ValueType>::value && has_unary_prod<ValueType>::value) &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function void y2_backward_subroutine(ValueType &y2_im1_, const ValueType &y2_i_, const ValueType &u_im1_)
        {
            // y2[i - 1] = y2[i - 1] * y2[i] + u[i - 1];
            auto y2_i = std::cbegin(y2_i_);
            auto u_im1 = std::cbegin(u_im1_);
            for (auto &&y2_im1 : y2_im1_) {
                (y2_im1 *= (*y2_i++)) += (*u_im1++);
            }
        }

        template<typename ValueType,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!(has_unary_plus<ValueType>::value && has_unary_prod<ValueType>::value) &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function void y2_backward_subroutine(ValueType &y2_im1_, const ValueType &y2_i_, const ValueType &u_im1_)
        {
            // y2[i - 1] = y2[i - 1] * y2[i] + u[i - 1];
            auto y2_i = std::cbegin(y2_i_);
            auto u_im1 = std::cbegin(u_im1_);
            for (auto &&y2_im1 : y2_im1_) {
                (y2_im1.second *= (y2_i->second)) += (u_im1->second);
                ++y2_i;
                ++u_im1;
            }
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_spline_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto slope_formula(const ValueType &v_iu, const ValueType &v_il, const Scalar &fac, const ValueType &y2_idx)
        {
            // slope = v[iu] - v[il] + fac * y2[idx], where fac = (+ | -) h * h * (1.0 / 6.0) and idx = (il | iu)
            return v_iu + fac * y2_idx - v_il;
                
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto slope_formula(ValueType v_iu, const ValueType &v_il_, const Scalar &fac, const ValueType &y2_idx_)
        {
            // slope = v[iu] - v[il] + fac * y2[idx], where fac = (+ | -) h * h * (1.0 / 6.0) and idx = (il | iu)
            auto y2_idx = std::cbegin(y2_idx_);
            auto v_il = std::cbegin(v_il_);
            for (auto &&s : v_iu) {
                s += (fac * (*y2_idx++) - *v_il++);
            }
            return v_iu;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto slope_formula(ValueType v_iu, const ValueType &v_il_, const Scalar &fac, const ValueType &y2_idx_)
        {
            // slope = v[iu] - v[il] + fac * y2[idx], where fac = (+ | -) h * h * (1.0 / 6.0) and idx = (il | iu)
            auto y2_idx = std::cbegin(y2_idx_);
            auto v_il = std::cbegin(v_il_);
            for (auto &&s : v_iu) {
                s.second += fac * (y2_idx->second) - v_il->second;
                ++y2_idx;
                ++v_il;
            }
            return v_iu;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 lookup::extrapolation ExtrapolationOption_ = extrapolation_option,
                 typename std::enable_if_t<has_spline_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::spline &&
                                           ExtrapolationOption_ == lookup::extrapolation::linear>* = nullptr>
        static constexpr pure_function auto yy_linextrap_formula(const ValueType &v_idx, const Scalar &fac, const ValueType &slope)
        {
            // yy[i] = v[idx] + fac * slope, where fac = (p - 1 | p) and idx = (iu | il)
            return v_idx + fac * slope;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 lookup::extrapolation ExtrapolationOption_ = extrapolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline &&
                                           ExtrapolationOption_ == lookup::extrapolation::linear>* = nullptr>
        static constexpr pure_function auto& yy_linextrap_formula(const ValueType &v_idx_, const Scalar &fac, ValueType &slope)
        {
            // yy[i] = v[idx] + fac * slope, where fac = (p - 1 | p) and idx = (iu | il)
            auto v_idx = std::cbegin(v_idx_);
            for (auto &&s : slope) {
                (s *= fac) += *v_idx++;
            }
            return slope;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 lookup::extrapolation ExtrapolationOption_ = extrapolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline &&
                                           ExtrapolationOption_ == lookup::extrapolation::linear>* = nullptr>
        static constexpr pure_function auto& yy_linextrap_formula(const ValueType &v_idx_, const Scalar &fac, ValueType &slope)
        {
            // yy[i] = v[idx] + fac * slope, where fac = (p - 1 | p) and idx = (iu | il)
            auto v_idx = std::cbegin(v_idx_);
            for (auto &&s : slope) {
                (s.second *= fac) += v_idx->second;
                ++v_idx;
            }
            return slope;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<has_spline_interp_ops<ValueType, Scalar>::value &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto yy_formula(const ValueType &v_il, const Scalar &p, const ValueType &v_iu,
                                                       const Scalar &smsq, const ValueType &y2_il, const Scalar &pmsq,
                                                       const ValueType &y2_iu, const Scalar &h2_div6)
        {
            // yy[i] = v[il] + p * (v[iu] - v[il]) + (smsq * y2[il] + pmsq * y2[iu]) * h * h * (1.0 / 6.0);
            return v_il + p * (v_iu - v_il) + h2_div6 * (smsq * y2_il + pmsq * y2_iu);
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           !(is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto yy_formula(ValueType v_il, const Scalar &p, const ValueType &v_iu_,
                                                       const Scalar &smsq, const ValueType &y2_il_, const Scalar &pmsq,
                                                       const ValueType &y2_iu_, const Scalar &h2_div6)
        {
            // yy[i] = v[il] + p * (v[iu] - v[il]) + (smsq * y2[il] + pmsq * y2[iu]) * h * h * (1.0 / 6.0);
            auto v_iu = std::cbegin(v_iu_);
            auto y2_il = std::cbegin(y2_il_);
            auto y2_iu = std::cbegin(y2_iu_);
            for (auto &&v : v_il) {
                // with +=, rounding errors make it return different results than the scalar version
                v = v + p * (*v_iu++ - v) + h2_div6 * (smsq * (*y2_il++) + pmsq * (*y2_iu++));
            }
            return v_il;
        }

        template<typename ValueType,
                 typename Scalar,
                 lookup::interpolation InterpolationOption_ = interpolation_option,
                 typename std::enable_if_t<!has_spline_interp_ops<ValueType, Scalar>::value &&
                                           (is_specialization<ValueType, std::map>::value || is_specialization<ValueType, std::unordered_map>::value) &&
                                           InterpolationOption_ == lookup::interpolation::spline>* = nullptr>
        static constexpr pure_function auto yy_formula(ValueType v_il, const Scalar &p, const ValueType &v_iu_,
                                                       const Scalar &smsq, const ValueType &y2_il_, const Scalar &pmsq,
                                                       const ValueType &y2_iu_, const Scalar &h2_div6)
        {
            // yy[i] = v[il] + p * (v[iu] - v[il]) + (smsq * y2[il] + pmsq * y2[iu]) * h * h * (1.0 / 6.0);
            auto v_iu = std::cbegin(v_iu_);
            auto y2_il = std::cbegin(y2_il_);
            auto y2_iu = std::cbegin(y2_iu_);
            for (auto &&v : v_il) {
                // with +=, rounding errors make it return different results than the scalar version
                v.second = v.second + p * (v_iu->second - v.second) + h2_div6 * (smsq * y2_il->second + pmsq * y2_iu->second);
                ++v_iu;
                ++y2_il;
                ++y2_iu;

            }
            return v_il;

        }


        template<lookup::interpolation InterpolationOption_, typename = void>
        struct extra_data_interpolation
        {
            template<typename... Args>
            static constexpr void build(Args&&...) {}

            static constexpr void clear() {}

        };

        template<lookup::interpolation InterpolationOption_>
        struct extra_data_interpolation<InterpolationOption_,
                                        std::void_t<std::enable_if_t<InterpolationOption_ == lookup::interpolation::spline> > >
        {
            template<typename Table, typename PrelookupType>
            constexpr void build(const Table &table,
                                 const std::array<std::size_t, rank> &grid_vector_sizes,
                                 const PrelookupType &first_prelookup)
            {
                // elts: the product of all following dimensions
                std::exclusive_scan(std::crbegin(grid_vector_sizes), std::crend(grid_vector_sizes),
                                    std::rbegin(elts), index_type{ 1 },
                                    std::multiplies<>{});

                // resize value buffers
                yyA_work.resize(elts.front());
                yyB_work.resize(elts.front());
                yy2_work.resize(elts.front());
                u_work.resize(*std::max_element(std::cbegin(grid_vector_sizes),
                                                std::cend(grid_vector_sizes)) - 1u);

                // establish the zero value
                auto t = std::cbegin(table);
                zero = zero_value<value_type, typename PrelookupType::grid_vector_type::value_type>(*t);

                // pre-calculate 2nd derivatives for the 1st dimension
                const auto len_0 = grid_vector_sizes.front();
                const auto max_index_0 = len_0 - 1u;
                y20_container.resize(len_0 * elts.front());
                auto y2 = std::begin(y20_container);
                const auto &first_grid_vector = first_prelookup.grid_vector();
                for (; t != std::cend(table); y2 += len_0, t += len_0) {
                    Lookup_table::cubic_spline_second_derivatives(first_grid_vector, max_index_0, t, y2, u_work, zero);
                }

            }

            constexpr void clear()
            {
                zero = value_type{};
                std::fill(std::begin(elts), std::end(elts), index_type{ 0 });
                y20_container.clear();
                yyA_work.clear();
                yyB_work.clear();
                yy2_work.clear();
                u_work.clear();

            }

            value_type zero;
            std::array<index_type, rank> elts;
            std::vector<value_type> y20_container;
            std::vector<value_type> yyA_work;
            std::vector<value_type> yyB_work;
            std::vector<value_type> yy2_work;
            std::vector<value_type> u_work;

        private:

            template<typename ValueType, typename Scalar,
                     typename std::enable_if_t<Lookup_table::has_spline_interp_ops<ValueType, Scalar>::value>* = nullptr>
            static constexpr pure_function ValueType zero_value(const ValueType &)
            {
                return ValueType{ 0 };
            }

            template<typename ValueType, typename Scalar,
                     typename std::enable_if_t<!Lookup_table::has_spline_interp_ops<ValueType, Scalar>::value &&
                                               !(is_specialization<ValueType, std::map>::value ||
                                                 is_specialization<ValueType, std::unordered_map>::value)>* = nullptr>
            static constexpr pure_function ValueType zero_value(ValueType token)
            {
                for (auto &&v : token) {
                    v *= Scalar{ 0 };
                }

                return token;

            }

            template<typename ValueType, typename Scalar,
                     typename std::enable_if_t<!Lookup_table::has_spline_interp_ops<ValueType, Scalar>::value &&
                                               (is_specialization<ValueType, std::map>::value ||
                                                is_specialization<ValueType, std::unordered_map>::value)>* = nullptr>
            static constexpr pure_function ValueType zero_value(ValueType token)
            {
                for (auto &&v : token) {
                    v.second *= Scalar{ 0 };
                }

                return token;

            }

        };

        mutable extra_data_interpolation<interpolation_option> _extra_data;

    };


    template<std::size_t Rank = 1u,
             typename PrelookupType = Prelookup<>,
             lookup::interpolation InterpolationOption = lookup::interpolation::linear,
             lookup::extrapolation ExtrapolationOption = lookup::extrapolation::allow,
             typename TableType = std::vector<double> >
    struct Homogeneous_lookup_table : Lookup_table<std::array<PrelookupType, Rank>, InterpolationOption, ExtrapolationOption, TableType>
    {
        //
        // Represents a lookup table (with "Rank" dimensions) that uses
        // the same prelookup type for every dimension.
        //

        using base_type = Lookup_table<std::array<PrelookupType, Rank>,
                                       InterpolationOption, ExtrapolationOption, TableType>;


        constexpr Homogeneous_lookup_table() = default;

        template<typename... Args>
        Homogeneous_lookup_table(Args&&... args) : base_type{ std::forward<Args>(args)... } {}

    };


    template<std::size_t Rank = 1u,
             lookup::interpolation InterpolationOption = lookup::interpolation::linear,
             lookup::extrapolation ExtrapolationOption = lookup::extrapolation::allow,
             typename TableType = std::vector<double> >
    struct Direct_lookup_table : Lookup_table<std::array<std::size_t, Rank>, InterpolationOption, ExtrapolationOption, TableType>
    {
        //
        // Represents a lookup table (with "Rank" dimensions) that
        // has no member prelookups, i.e., it will always be queried
        // with "prelookup_result" inputs, which may be obtained
        // from external "tr7::Prelookup" objects. In these tables,
        // the "_prelookups" data member will only be an std::array
        // holding the size of the table along each dimension.
        //

        using base_type = Lookup_table<std::array<std::size_t, Rank>,
                                       InterpolationOption, ExtrapolationOption, TableType>;


        static_assert(InterpolationOption != lookup::interpolation::spline,
                      "tr7::Direct_lookup_table: the \"spline\" interpolation option"
                      " is not supported for this type of lookup table.");


        constexpr Direct_lookup_table() = default;

        template<typename... Args>
        Direct_lookup_table(Args&&... args) : base_type{ std::forward<Args>(args)... } {}

    };


}


#endif
