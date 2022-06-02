#ifndef F16_NGUYEN_F16_NGUYEN_ENGINE_H
#define F16_NGUYEN_F16_NGUYEN_ENGINE_H
#pragma once


#include "tr7/lookup.h"

#include "F16_Nguyen/modelglue.h"
#include "F16_Nguyen/contiguous_ranges.h"
#include "F16_Nguyen/read_csv_table.h"


//
// Defines the "F16_Nguyen_engine" class, which implements
// the engine model explained in reference "Simulator study
// of stall/post-stall characteristics of a fighter airplane with
// relaxed longitudinal static stability", L.T.Nguyen et al.,
// NASA-TP-1538, 1979.
//


namespace F16_Nguyen
{
    class F16_Nguyen_engine
    {
        //
        // A class to evaluate Nguyen's F16 engine model.
        // It has one state variable ("P3_percent") and
        // several outputs, with the main one being the
        // engine's thrust.
        //

    public:

        using data_type = double;


        struct default_parameters
        {
            static constexpr data_type engine_angular_momentum_kgm2ps = 216.9;
        };


        DECLARE_MODEL_INPUTS(throttle_percent, Mach, ZP_m)

        DECLARE_MODEL_STATES(P3_percent)

        DECLARE_MODEL_SIGNALS(dataset_output, P1_percent,
                                              Tidle_N, Tmil_N, Tmax_N)

        DECLARE_MODEL_OUTPUTS_DERIVED_FROM(dataset_output, engine_thrust_N,
                                                           engine_angular_momentum_kgm2ps,
                                                           P2_percent,
                                                           invtau_1ps)


        F16_Nguyen_engine() = default;

        template<typename... Args>
        F16_Nguyen_engine(Args&&... args)
        {
            build(std::forward<Args>(args)...);
        }

        void build(const std::string &path_dataset);


        template<typename StatesType, typename T>
        constexpr outputs_type<T> outputs(const StatesType &x,
                                          const T &throttle_percent,
                                          const T &Mach,
                                          const T &ZP_m = T{ 0 }) const;

        template<typename StatesType, typename InputsType>
        constexpr auto outputs(const StatesType &x,
                               const InputsType &u) const
        {
            ensure_minimum_size<_num_inputs>(u);
            return  outputs(x,
                            u[input_names::throttle_percent],
                            u[input_names::Mach],
                            u[input_names::ZP_m]);

        }


        template<typename StatesType, typename OutputsType,
                 typename T = tr7::typeof_bracket_operator<StatesType> >
        static constexpr states_type<T> derivatives(const StatesType &x,
                                                    const OutputsType &y);

        template<typename StatesType,
                 typename T = tr7::typeof_bracket_operator<StatesType> >
        static constexpr states_type<T> state_limiters(const StatesType &x);


        template<typename T>
        constexpr dataset_outputs_type<T> dataset(const T &throttle_percent,
                                                  const T &Mach,
                                                  const T &ZP_m) const;

        template<typename InputsType>
        constexpr auto dataset(const InputsType &u) const
        {
            ensure_minimum_size<_num_inputs>(u);
            return  dataset(u[input_names::throttle_percent],
                            u[input_names::Mach],
                            u[input_names::ZP_m]);

        }


        template<typename T>
        constexpr T invtau_1ps(const T &P2minusP3_percent) const;


    private:

        template<typename T>
        constexpr dataset_outputs_type<T> lookup_dataset(const T &throttle_percent,
                                                         const T &Mach,
                                                         const T &ZP_m) const;


        struct luts
        {
            static constexpr auto prelookup_option = tr7::lookup::prelookup::linsearch;
            static constexpr auto use_previous_index = true;

            using prelookup_type = tr7::Prelookup<prelookup_option,
                                                  use_previous_index,
                                                  std::vector<data_type> >;

            static constexpr auto interpolation_option = tr7::lookup::interpolation::linear;
            static constexpr auto extrapolation_option = tr7::lookup::extrapolation::nearest;

            using lut1d_type = tr7::Homogeneous_lookup_table<1u,
                                                             prelookup_type,
                                                             interpolation_option,
                                                             extrapolation_option,
                                                             std::vector<data_type> >;

            using lut2d_type = tr7::Homogeneous_lookup_table<2u,
                                                             prelookup_type,
                                                             interpolation_option,
                                                             extrapolation_option,
                                                             std::vector<std::array<data_type, 3u> > >;

            // plain 1D luts
            lut1d_type P1_percent; // "f(throttle)"
            lut1d_type invtau_1ps; // "f(P2 - P3)"

            // 2D lut, which is "f(Mach, ZP_m)"
            lut2d_type Tidle_Tmil_Tmax_N;

        } _luts;

    };


    //
    // "F16_Nguyen_engine" impl.
    //

    inline void F16_Nguyen_engine::build(const std::string &path_dataset)
    {
        //
        // Reads the "NGUYEN F16 ENGINE-DATASET" "*.csv" files living in
        // "path_dataset", and sets up the lookup tables that we'll use
        // to evaluate the aircraft's engine model.
        //

        // helper function to read a plain 1D lut
        const auto read_lut1d = [&path_dataset](const auto &lutname,
                                                const auto &breaks_filename, const auto &table_filename)
        {
            auto breaks = std::get<0>(read_csv_table(path_dataset + "/" + breaks_filename));
            auto [t, num_rows, num_columns] = read_csv_table(path_dataset + "/" + table_filename);

            if (num_rows != breaks.size() || num_columns != 1u) {
                TR7_THROW_RUNTIME_ERROR(std::string{ "F16_Nguyen::F16_Nguyen_engine::build: inconsistent sizes while"
                                                     " reading \"" } + lutname + "\".")

            }

            return luts::lut1d_type(std::move(t), std::move(breaks));

        };


        // read the object's 1D luts
        _luts.P1_percent = read_lut1d("luts::P1_percent", "throttle_breaks_percent.csv", "P1_throttle.csv");
        _luts.invtau_1ps = read_lut1d("luts::invtau_1ps", "P2minusP3_breaks_percent.csv", "invtau_P2minusP3.csv");


        // read the "Tidle_Tmin_Tmax_N" table
        {
            auto Mach_breaks = std::get<0>(read_csv_table(path_dataset + "/Mach_breaks.csv"));
            auto ZP_breaks_m = std::get<0>(read_csv_table(path_dataset + "/alt_breaks_m.csv"));

            const auto size_x = Mach_breaks.size();
            const auto size_y = ZP_breaks_m.size();

            const auto [Tidle_N, num_rows_Tidle_N, num_columns_Tidle_N] = read_csv_table(path_dataset + "/Tidle_Mach_alt.csv");
            const auto [Tmil_N, num_rows_Tmil_N, num_columns_Tmil_N] = read_csv_table(path_dataset + "/Tmil_Mach_alt.csv");
            const auto [Tmax_N, num_rows_Tmax_N, num_columns_Tmax_N] = read_csv_table(path_dataset + "/Tmax_Mach_alt.csv");

            if (num_rows_Tidle_N != size_x || num_rows_Tmil_N != size_x || num_rows_Tmax_N != size_x ||
                num_columns_Tidle_N != size_y || num_columns_Tmil_N != size_y || num_columns_Tmax_N != size_y) {

                TR7_THROW_RUNTIME_ERROR("F16_Nguyen::F16_Nguyen_engine::build: inconsistent sizes while"
                                        " reading \"luts::Tidle_Tmil_Tmax_N\".")

            }

            typename luts::lut2d_type::table_type t(size_x * size_y);
            for (auto count = 0u; count < size_x * size_y; ++count) {
                t[count] = { Tidle_N[count], Tmil_N[count], Tmax_N[count] };
            }

            _luts.Tidle_Tmil_Tmax_N.build(std::move(t), std::move(Mach_breaks), std::move(ZP_breaks_m));

        }

    }


    template<typename StatesType, typename T>
    constexpr auto F16_Nguyen_engine::outputs(const StatesType &x,
                                              const T &throttle_percent,
                                              const T &Mach,
                                              const T &ZP_m) const -> outputs_type<T>
    {
        //
        // Evaluates the model's output equation, i.e.,
        // "y = G(x, u, ...)", and returns the model's "y".
        //

        ensure_minimum_size<_num_states>(x);


        // calculate the regime
        const auto P3_nondim = x[state_names::P3_percent] / T{ 50 };
        const auto high_regime = P3_nondim >= T{ 1 };


        // evaluate the model's dataset (P1 and regime-wise thrusts)
        const auto [P1_percent,
                    Tidle_N, Tmil_N, Tmax_N] = dataset(throttle_percent,
                                                       Mach,
                                                       ZP_m);


        // convert (P1, P3) into (P2, invtau)
        T P2_percent;
        T invtau_1ps;
        if (P1_percent >= T{ 50 }) {
            if (high_regime) {
                P2_percent = P1_percent;
                invtau_1ps = T{ 5 };

            }
            else {
                P2_percent = T{ 60 };
                invtau_1ps = this->invtau_1ps(P2_percent - x[state_names::P3_percent]);

            }

        }
        else if (high_regime) {
            P2_percent = T{ 40 };
            invtau_1ps = T{ 5 };

        }
        else {
            P2_percent = P1_percent;
            invtau_1ps = this->invtau_1ps(P2_percent - x[state_names::P3_percent]);

        }


        // calculate the final engine thrust
        const auto engine_thrust_N = high_regime ? Tmil_N + (Tmax_N - Tmil_N) * (P3_nondim - T{ 1 }) :
                                                   Tidle_N + (Tmil_N - Tidle_N) * P3_nondim;


        // the engine's angular momentum is taken constant in the reference (?)
        const auto engine_angular_momentum_kgm2ps = T{ default_parameters::engine_angular_momentum_kgm2ps };


        // return the model's outputs
        return outputs_type<T>{ P1_percent,
                                Tidle_N, Tmil_N, Tmax_N,
                                engine_thrust_N,
                                engine_angular_momentum_kgm2ps,
                                P2_percent,
                                invtau_1ps };

    }


    template<typename StatesType, typename OutputsType,
             typename T>
    constexpr auto F16_Nguyen_engine::derivatives(const StatesType &x,
                                                  const OutputsType &y) -> states_type<T>
    {
        //
        // Evaluates the model's state equation, i.e.,
        // "xdot = F(x, y, u, ...)", and returns the
        // model's "xdot".
        //

        ensure_minimum_size<_num_states>(x);
        ensure_minimum_size<_num_outputs>(y);


        // calculate P3dot and return it
        return states_type<T>{ y[output_names::invtau_1ps] * (y[output_names::P2_percent] - x[state_names::P3_percent]) };

    }


    template<typename StatesType,
             typename T>
    constexpr auto F16_Nguyen_engine::state_limiters(const StatesType &x) -> states_type<T>
    {
        //
        // Evaluates the model's state limiter equation, i.e.,
        // "x = SL(x, ...)", and returns the model's limited "x".
        //

        ensure_minimum_size<_num_states>(x);


        // P3 is a percentage in [0, 100].
        if (x[state_names::P3_percent] <= T{ 0 }) {
            return states_type<T>{ 0 };
        }
        else if (x[state_names::P3_percent] >= T{ 100 }) {
            return states_type<T>{ 100 };
        }
        else {
            return states_type<T>{ x[state_names::P3_percent] };
        }

    }



    template<typename T>
    constexpr auto F16_Nguyen_engine::dataset(const T &throttle_percent,
                                              const T &Mach,
                                              const T &ZP_m) const -> dataset_outputs_type<T>
    {
        //
        // Evaluates the dataset's coefficients, returning a
        // central difference approximation when we're in AD.
        //

        if constexpr (!tr7::is_AD_scalar<T>::value) {
            return lookup_dataset(throttle_percent, Mach, ZP_m);

        }
        else {
            // the AD with a linear lookup table behaves as a forward difference,
            // we need to add the influence of the backward one to get a central
            // result. We'll only do it for the thrust lut, not for the throttle
            // one (since that one is not physical, i.e. the engine system is
            // designed to have an asymmetrical response w.r.t. the throttle input)
            auto douts = lookup_dataset(throttle_percent, Mach, ZP_m);

            const auto Mach_numerical = tr7::AD_force_scalar_value(Mach);
            const auto ZP_numerical_m = tr7::AD_force_scalar_value(ZP_m);

            constexpr auto delta_numjac = decltype(Mach_numerical){ 1e-9 };
            const auto Tidle_Tmil_Tmax_Mach_minus = _luts.Tidle_Tmil_Tmax_N(Mach_numerical - delta_numjac, ZP_numerical_m);
            const auto Tidle_Tmil_Tmax_ZP_minus = _luts.Tidle_Tmil_Tmax_N(Mach_numerical, ZP_numerical_m - delta_numjac);

            const auto DMach_div_delta = (Mach - Mach_numerical) / delta_numjac;
            const auto DZP_div_delta = (ZP_m - ZP_numerical_m) / delta_numjac;

            const auto Tidle_numerical_N = tr7::AD_force_scalar_value(douts[dataset_output_names::Tidle_N]);
            const auto DTidle_Mach = Tidle_numerical_N - Tidle_Tmil_Tmax_Mach_minus[0];
            const auto DTidle_ZP = Tidle_numerical_N - Tidle_Tmil_Tmax_ZP_minus[0];
            douts[dataset_output_names::Tidle_N] = T{ 0.5 } * douts[dataset_output_names::Tidle_N] +
                                                    decltype(Mach_numerical){ 0.5 } * (Tidle_numerical_N +
                                                                                       DTidle_Mach * DMach_div_delta +
                                                                                       DTidle_ZP * DZP_div_delta);

            const auto Tmil_numerical_N = tr7::AD_force_scalar_value(douts[dataset_output_names::Tmil_N]);
            const auto DTmil_Mach = Tmil_numerical_N - Tidle_Tmil_Tmax_Mach_minus[1];
            const auto DTmil_ZP = Tmil_numerical_N - Tidle_Tmil_Tmax_ZP_minus[1];
            douts[dataset_output_names::Tmil_N] = T{ 0.5 } * douts[dataset_output_names::Tmil_N] +
                                                    decltype(Mach_numerical){ 0.5 } * (Tmil_numerical_N +
                                                                                       DTmil_Mach * DMach_div_delta +
                                                                                       DTmil_ZP * DZP_div_delta);

            const auto Tmax_numerical_N = tr7::AD_force_scalar_value(douts[dataset_output_names::Tmax_N]);
            const auto DTmax_Mach = Tmax_numerical_N - Tidle_Tmil_Tmax_Mach_minus[2];
            const auto DTmax_ZP = Tmax_numerical_N - Tidle_Tmil_Tmax_ZP_minus[2];
            douts[dataset_output_names::Tmax_N] = T{ 0.5 } * douts[dataset_output_names::Tmax_N] +
                                                    decltype(Mach_numerical){ 0.5 } * (Tmax_numerical_N +
                                                                                       DTmax_Mach * DMach_div_delta +
                                                                                       DTmax_ZP * DZP_div_delta);

            return douts;

        }

    }


    template<typename T>
    constexpr auto F16_Nguyen_engine::lookup_dataset(const T &throttle_percent,
                                                     const T &Mach,
                                                     const T &ZP_m) const -> dataset_outputs_type<T>
    {
        //
        // Evaluates the model's lookup table in which we
        // store the regime-wise thrusts, and the
        // "throttle -> P1" curve.
        //

        const auto P1_percent = _luts.P1_percent(throttle_percent);
        const auto Tidle_Tmil_Tmax = _luts.Tidle_Tmil_Tmax_N(Mach, ZP_m);

        return dataset_outputs_type<T>{ T{ P1_percent },
                                        T{ std::get<0>(Tidle_Tmil_Tmax) },
                                        T{ std::get<1>(Tidle_Tmil_Tmax) },
                                        T{ std::get<2>(Tidle_Tmil_Tmax) } };

    }


    template<typename T>
    constexpr T F16_Nguyen_engine::invtau_1ps(const T &P2minusP3_percent) const
    {
        //
        // Evaluates the "(P2 - P3) -> 1 / tau_s" model curve.
        // BEWARE: the actual value of "tau_s" depends on the
        // engine's regime, this curve is only queried on some
        // cases.
        //

        return T{ _luts.invtau_1ps(P2minusP3_percent) };

    }


}


#endif