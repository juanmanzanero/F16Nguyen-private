#ifndef F16_NGUYEN_MODELGLUE
#define F16_NGUYEN_MODELGLUE
#pragma once


#include <array>
#include <algorithm>
#include <string>
#include <sstream>

#include "tr7/tr7tuple.h"


//
// Defines "detail" templates, functions & macros
// that we'll use to prepare & combine models, i.e.
// the enums they use to define their scalar
// outputs/states, the std arrays that will hold
// them, etc.
//


namespace F16_Nguyen
{
        //
        // Macros to declare a model's SCALAR states/outputs/inputs/inflags
        // (or other signals that don't intervene directly in the
        // state-space interface, but may be interesting) as enums,
        // also providing:
        //
        //    * A helper static size_t member "_num_(states/outputs/inputs/inflags)",
        //      equal to the number of elements, and a static function
        //      "num_(states/outputs/inputs/inflags)" that returns this value.
        //
        //    * A helper template typedef "(states/outputs/inputs/inflags)_type"
        //      to have them as an std::array of size "num_(states/outputs/inputs/inflags)",
        //      and value_type "T".
        //
        //    * A helper function "(output/state/input/inflag)_names2str" that
        //      returns their names as an std::array of std::strings.
        //
        // NOTE: we also provide versions that derived models should
        // use (THEY ONLY WORK FOR A SINGLE BASE).
        //
        // Example:
        //
        //     struct Basic3Outputs
        //     {
        //         DECLARE_MODEL_OUTPUTS(out0, out1, out2)
        //         ... // member functions, etc
        //     };
        //
        //     struct My_model
        //     {
        //         DECLARE_MODEL_STATES(x0, x1, x2, x3, x4) // _num_states = 4; states_type<T> = std::array<T, 4>; state_names2str() returns an std::array<std::string, 4>{ "x1", "x2", "x3", "x4" }
        //         DECLARE_MODEL_DERIVED_OUTPUTS(Basic3Outputs, my_y0, my_y1) // _num_outputs = 5; outputs_type<T> = std::array<T, 5>; output_names2str() returns an std::array<std::string, 5>{ "out0", "out1", "out2", "my_y0", "my_y1" }
        //         DECLARE_MODEL_INPUTS(u0, u1) // _num_inputs = 2; inputs_type<T> = std::array<T, 2>; input_names2str() returns an std::array<std::string, 2>{ "u0", "u1" }
        //         DECLARE_MODEL_INFLAGS(mode, submode) // ("inflags" are integers that make the model behave in some way)
        //         ... // member functions, etc
        //     };
        //

#define DETAIL_STRINGIFY_LINE(...) #__VA_ARGS__

#define DETAIL_DECLARE_INTERFACE_ENUM(whatever, \
                                      FIRSTNAME, ...) template<std::size_t Begin> \
                                                      struct whatever##_names_T \
                                                      { \
                                                          enum { FIRSTNAME = Begin, ##__VA_ARGS__, _count }; \
                                                      }; \
                                                      \
                                                      using whatever##_names = whatever##_names_T<0u>; \
                                                      static constexpr std::size_t _num_##whatever##s = whatever##_names::_count; \
                                                      static constexpr std::size_t num_##whatever##s() { return _num_##whatever##s; } \
                                                      \
                                                      template<typename T> \
                                                      using whatever##s_type = std::array<T, _num_##whatever##s>; \
                                                      \
                                                      static inline const std::array<std::string, _num_##whatever##s>& whatever##_names2str() \
                                                      { \
                                                          static const auto arr = detail::ss2strarray<_num_##whatever##s>(std::stringstream(DETAIL_STRINGIFY_LINE(FIRSTNAME, ##__VA_ARGS__))); \
                                                          return arr; \
                                                      }


#define DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(whatever, \
                                              BASETYPE, FIRSTNAME, ...) template<std::size_t Begin> \
                                                                        struct whatever##_names_T : BASETYPE##_names_T<Begin> \
                                                                        { \
                                                                            using base_names = BASETYPE##_names_T<Begin>; \
                                                                            enum { FIRSTNAME = base_names::_count, ##__VA_ARGS__, _count }; \
                                                                            static constexpr std::size_t _local_count = _count - base_names::_count; \
                                                                        }; \
                                                                        \
                                                                        using whatever##_names = whatever##_names_T<0u>; \
                                                                        static constexpr std::size_t _num_##whatever##s = whatever##_names::_count; \
                                                                        static constexpr std::size_t num_##whatever##s() { return _num_##whatever##s; } \
                                                                        \
                                                                        template<typename T> \
                                                                        using whatever##s_type = std::array<T, _num_##whatever##s>; \
                                                                        \
                                                                        static inline const std::array<std::string, _num_##whatever##s>& whatever##_names2str() \
                                                                        { \
                                                                            static const auto arr = detail::strarray_cat(BASETYPE##_names2str(), \
                                                                                                                         detail::ss2strarray<whatever##_names::_local_count>(std::stringstream(DETAIL_STRINGIFY_LINE(FIRSTNAME, ##__VA_ARGS__)))); \
                                                                            return arr; \
                                                                        }


// declare an enum of states/outputs/inputs/inflags/signals
#define DECLARE_MODEL_STATES(FIRSTNAME, ...) DETAIL_DECLARE_INTERFACE_ENUM(state, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_OUTPUTS(FIRSTNAME, ...) DETAIL_DECLARE_INTERFACE_ENUM(output, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_INPUTS(FIRSTNAME, ...) DETAIL_DECLARE_INTERFACE_ENUM(input, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_INFLAGS(FIRSTNAME, ...) DETAIL_DECLARE_INTERFACE_ENUM(inflag, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_SIGNALS(ENUMNAME_IN_SINGULAR, FIRSTNAME, ...) DETAIL_DECLARE_INTERFACE_ENUM(ENUMNAME_IN_SINGULAR, FIRSTNAME, ##__VA_ARGS__)

// declare an enum of states/outputs/inputs/inflags/signals that inherits
// from an enum of the same name declared in another class
#define DECLARE_MODEL_DERIVED_STATES(CLASSNAME, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(state, CLASSNAME::state, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_DERIVED_OUTPUTS(CLASSNAME, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(output, CLASSNAME::output, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_DERIVED_INPUTS(CLASSNAME, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(input, CLASSNAME::input, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_DERIVED_INFLAGS(CLASSNAME, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(inflag, CLASSNAME::inflag, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_DERIVED_SIGNALS(ENUMNAME_IN_SINGULAR, CLASSNAME, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(ENUMNAME_IN_SINGULAR, CLASSNAME::ENUMNAME_IN_SINGULAR, FIRSTNAME, ##__VA_ARGS__)

// declare a list of states/outputs/inputs/inflags/signals that inherits
// from another enum with the given name
#define DECLARE_MODEL_STATES_DERIVED_FROM(BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(state, BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_OUTPUTS_DERIVED_FROM(BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(output, BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_INPUTS_DERIVED_FROM(BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(input, BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_INFLAGS_DERIVED_FROM(BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(inflag, BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ##__VA_ARGS__)
#define DECLARE_MODEL_SIGNALS_DERIVED_FROM(THIS_ENUMNAME_IN_SINGULAR, BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ...) DETAIL_DECLARE_DERIVED_INTERFACE_ENUM(THIS_ENUMNAME_IN_SINGULAR, BASE_ENUMNAME_IN_SINGULAR, FIRSTNAME, ##__VA_ARGS__)


    namespace detail
    {
        //
        // Convert states/outputs/inputs of a model into std::strings (and
        // concatenate those of several submodels). I have to do this at
        // runtime, and using std::string instead of const char*, because
        // I've not been able to create a recursive "STRINGIFY" macro that
        // acts on a comma-separated list of names... maybe some day...
        //

        template<std::size_t N>
        inline std::array<std::string, N> ss2strarray(std::stringstream &&ss)
        {
            //
            // Splits an std::stringstream holding a comma-separated
            // list of labels (comma + space, e.g. "hola, kpasa, tio")
            // into an "std::array<std::string, N>", where "N" should
            // be smaller or equal to the number of labels in "ss".
            //

            std::array<std::string, N> ret;

            auto r = std::begin(ret);
            while (ss.good() && r != std::end(ret)) {
                std::string s;
                std::getline(ss, s, ',');
                *r++ = std::move(s);
                ss.ignore();

            }

            if (std::any_of(std::cbegin(ret), std::cend(ret), [](auto &name) { return name.empty(); })) {
                TR7_THROW_RUNTIME_ERROR("F6_Nguyen::detail::ss2strarray: input labels could not be correctly"
                                        " converted into an array of strings.")

            }

            return ret;

        }

        template<typename Arr, typename... Arrs>
        struct sum_of_array_sizes
        {
            static constexpr std::size_t value = std::tuple_size_v<std::decay_t<Arr>> + sum_of_array_sizes<Arrs...>::value;
        };

        template<typename Arr>
        struct sum_of_array_sizes<Arr>
        {
            static constexpr std::size_t value = std::tuple_size_v<std::decay_t<Arr>>;
        };

        template<typename... Arrs>
        inline constexpr std::size_t sum_of_array_sizes_v = sum_of_array_sizes<Arrs...>::value;

        template<typename... Arrs>
        inline std::array<std::string, sum_of_array_sizes_v<Arrs...>> strarray_cat(Arrs&&... arrs)
        {
            //
            // Concatenates a pack of "std::array<std::string, Ni>"
            // into a single "std::array<std::string, sum(Ni)>"
            // holding all of their values.
            //

            std::array<std::string, sum_of_array_sizes_v<Arrs...>> ret;
            auto r = std::begin(ret);

            tr7::tuple_for(std::forward_as_tuple(arrs...), [&r](const auto &arr)
                                                           {
                                                               std::copy(std::cbegin(arr), std::cend(arr), r);
                                                               r += std::tuple_size_v<std::decay_t<decltype(arr)>>;
                                                           });

            return ret;

        }


        //
        // Concatenate the state/output/input names of a
        // collection of submodels, in a submodel-wise
        // contiguous fashion.
        //

        template<template<std::size_t, typename> class GetNames_t, std::size_t Begin, typename Model, typename... RestOfModels>
        struct glue_names_recursive : GetNames_t<Begin, Model>, glue_names_recursive<GetNames_t, GetNames_t<Begin, Model>::_count, RestOfModels...>
        {
            using _model = Model;
            using _names = GetNames_t<Begin, Model>;

            static constexpr std::size_t _begin = Begin;
            static constexpr std::size_t _end = _names::_count;

            using _next = glue_names_recursive<GetNames_t, _end, RestOfModels...>;

        };

        template<template<std::size_t, typename> class GetNames_t, std::size_t Begin, typename Model>
        struct glue_names_recursive<GetNames_t, Begin, Model> : GetNames_t<Begin, Model>
        {
            using _model = Model;
            using _names = GetNames_t<Begin, Model>;

            static constexpr std::size_t _begin = Begin;
            static constexpr std::size_t _end = _names::_count;

        protected:

            // (_total_begin is guaranteed to be 0u, and _total_end == _end at this level)
            static constexpr std::size_t _num_names = _end;

        };

        template<template<std::size_t, typename> class GetNames_t, typename... Models>
        struct glue_names : glue_names_recursive<GetNames_t, 0u, Models...>
        {
            using _root = glue_names_recursive<GetNames_t, 0u, Models...>;
            using _root::_num_names;

        private:

            using _root::_model;
            using _root::_begin;
            using _root::_end;

        };

        template<std::size_t Begin, class Model>
        using get_state_names_t = typename Model::template state_names_T<Begin>;

        template<std::size_t Begin, class Model>
        using get_output_names_t = typename Model::template output_names_T<Begin>;

        template<std::size_t Begin, class Model>
        using get_inflag_names_t = typename Model::template inflag_names_T<Begin>;

        template<typename... Models>
        using glue_output_names = glue_names<get_output_names_t, Models...>;

        template<typename... Models>
        using glue_state_names = glue_names<get_state_names_t, Models...>;

        template<typename... Models>
        using glue_inflag_names = glue_names<get_inflag_names_t, Models...>;

    }


}


#endif