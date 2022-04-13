#ifndef F16_NGUYEN_F16_NGUYEN_TRIM_H
#define F16_NGUYEN_F16_NGUYEN_TRIM_H
#pragma once


#include "F16_Nguyen/anemometry.h"
#include "F16_Nguyen/F16_Nguyen_aerodynamics.h"
#include "F16_Nguyen/F16_Nguyen_engine.h"
#include "F16_Nguyen/F16_Nguyen_flat_Earth_body_dynamics.h"


//
// Defines the "F16_Nguyen_trim" class, that takes
// an "F16_Nguyen_plant" and trims it.
//


namespace F16_Nguyen
{
    class F16_Nguyen_trim
    {
        //
        // Provides trim points for the
        // "F16_Nguyen_plant".
        //

    public:

        DECLARE_MODEL_INFLAGS(steady_trim, straight_and_level_trim, disable_lef)

        DECLARE_MODEL_INPUTS(KCAS, ZP_ft, DISA_K, flight_path_angle_deg, mass_kg, xcg_per_MAC)

        DECLARE_MODEL_OUTPUTS(cost_Fx, cost_Fy, cost_Fz, cost_Mx, cost_My, cost_Mz,
                              pitch_deg, dh_deg, aoa_deg, dlef_deg, P3_percent, throttle_percent)


        template<typename T = int>
        static constexpr inflags_type<T> default_inflags();

        template<typename T = double>
        static constexpr inputs_type<T> default_inputs();


        template<typename F16NguyenPlantType,
                 typename InflagsType, typename InputsType,
                 typename T = tr7::typeof_bracket_operator<InputsType> >
        static std::pair<bool, outputs_type<T> > trim_plant(const F16NguyenPlantType &plant,
                                                            const InflagsType &inflags,
                                                            const InputsType &u)
        {
            //
            // Trims the plant attending to the given inputs & inflags.
            //

            ensure_minimum_size<_num_inflags>(inflags);
            ensure_minimum_size<_num_inputs>(u);

            if (static_cast<int>(inflags[inflag_names::steady_trim]) &&
                static_cast<int>(inflags[inflag_names::straight_and_level_trim])) {

                return straight_and_level_trim(plant, inflags, u);

            }
            else {
                TR7_THROW_RUNTIME_ERROR("F16_Nguyen_trim::trim_plant: unsupported trim mode.")

            }

        }


    private:

        template<typename F16NguyenPlantType,
                 typename InflagsType, typename InputsType,
                 typename T = tr7::typeof_bracket_operator<InputsType> >
        static std::pair<bool, outputs_type<T> > straight_and_level_trim(const F16NguyenPlantType &plant,
                                                                         const InflagsType &inflags,
                                                                         const InputsType &u);


        template<typename F16NguyenPlantType, typename T>
        static std::tuple<bool, T, T, T> engine_trim(const F16NguyenPlantType &plant,
                                                     const T &necessary_thrust_N,
                                                     const T &Mach, const T &ZP_m);


        template<typename InflagsType, typename T>
        static constexpr T stationary_lef_schedule(const T &aoa_deg, const T &q_div_p,
                                                   const InflagsType &inflags);

        template<typename F16NguyenPlantType,
                 typename InflagsType, typename T>
        static constexpr T aoa_wind_point_seed(const F16NguyenPlantType &plant,
                                               const T &CW,
                                               const T &q_div_p,
                                               const F16_Nguyen_aerodynamics::force_and_moment_inputs_type<T> &basic_aeroforce_and_moment_inputs,
                                               const InflagsType &inflags,
                                               const T &Daoa_numjac_deg = T{ 1e-9 },
                                               std::size_t max_iters_point = 2u);

        template<typename F16NguyenPlantType, typename T>
        static constexpr T stabilator_point_seed(const F16NguyenPlantType &plant,
                                                 const F16_Nguyen_aerodynamics::force_and_moment_inputs_type<T> &basic_aeroforce_and_moment_inputs,
                                                 const T &Ddh_numjac_deg = T{ 1e-9 },
                                                 std::size_t max_iters_point = 2u);

    };


    //
    // "F16_Nguyen_trim" impl.
    //

    template<typename T>
    constexpr auto F16_Nguyen_trim::default_inflags() -> inflags_type<T>
    {
        //
        // Returns an array of default inflags so that the user
        // doesn't need to set all of them every time they
        // call the object.
        //

        constexpr auto steady_trim = T{ 1 };
        constexpr auto straight_and_level_trim = T{ 1 };
        constexpr auto disable_lef = T{ 0 };

        return inflags_type<T>{ steady_trim, straight_and_level_trim, disable_lef };

    }


    template<typename T>
    constexpr auto F16_Nguyen_trim::default_inputs() -> inputs_type<T>
    {
        //
        // Returns an array of default inputs so that the user
        // doesn't need to set all of them every time they
        // call the object.
        //

        constexpr auto KCAS = T{ 200 };
        constexpr auto ZP_ft = T{ 0 };
        constexpr auto DISA_K = T{ 0 };
        constexpr auto flight_path_angle_deg = T{ 0 };
        constexpr auto mass_kg = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::mass_kg };
        constexpr auto xcg_per_MAC = T{ F16_Nguyen_aerodynamics::default_parameters::reference_xcg_per_MAC };

        return inputs_type<T>{ KCAS, ZP_ft, DISA_K, flight_path_angle_deg, mass_kg, xcg_per_MAC };

    }


    template<typename F16NguyenPlantType,
             typename InflagsType, typename InputsType,
             typename T>
    inline auto F16_Nguyen_trim::straight_and_level_trim(const F16NguyenPlantType &plant,
                                                         const InflagsType &inflags,
                                                         const InputsType &u) -> std::pair<bool, outputs_type<T> >
    {
        //
        // Calculates the trimmed flight pitch of the "F16_Nguyen_plant", by solving
        // the simplified longitudinal balance that assumes only two control
        // surfaces (horizontal stabilator and leading edge flap):
        //
        //   q * S * CZ(aoa, dh, dlef) / W + cos(theta)      = 0,
        //   CM(aoa, dh, dlef)                               = 0,
        //   aoa + gamma - theta                             = 0,
        //   dlef_deg - 1.38 * aoa_deg + 9.05 * q / p - 1.45 = 0, (stationary leading edge flap scheduling law, from the "Nguyen" reference)
        //
        // After balancing the aerodynamics, we can extract the engine's state
        // "P3" and input "throttle" by solving:
        //
        //   thrust(P3) / W - sin(theta) + q * S * CX(aoa, dh, dlef) / W = 0,
        //   P3dot(throttle)                                             = 0.
        //
        // NOTE: this trim only evaluates the plant's aero-force and
        // moment coefficients, not its state-space functions, and considers
        // an engine totally aligned with the x-body-axis.
        //

        // basic S & L trim constants
        const auto ZP_m = tr7::ft2m(u[input_names::ZP_ft]);
        const auto atmos = isa_atmosphere(ZP_m, u[input_names::DISA_K]);
        const auto TAS_mps = cas2tas(tr7::kn2mps(u[input_names::KCAS]), ZP_m, u[input_names::DISA_K]);
        const auto dynamic_pressure_Pa = 0.5 * atmos.air_density_kgpm3 * TAS_mps * TAS_mps;
        const auto q_div_p = dynamic_pressure_Pa / atmos.air_pressure_Pa;
        const auto weight_N = T{ g0_mps2 } * u[input_names::mass_kg];
        const auto qS_div_W = dynamic_pressure_Pa * T{ F16_Nguyen_aerodynamics::default_parameters::wing_surface_m2 } / weight_N;


        // create the basic inputs to evaluate the aerodynamic model
        using aeroforce_and_moment_inputs_type = F16_Nguyen_aerodynamics::force_and_moment_inputs_type<T>;
        using aeroforce_and_moment_input_names = F16_Nguyen_aerodynamics::force_and_moment_input_names;

        aeroforce_and_moment_inputs_type basic_aeroforce_and_moment_inputs{};
        basic_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::TAS_mps] = TAS_mps;
        basic_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::xcg_per_MAC] = u[input_names::xcg_per_MAC];
        basic_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::MAC_m] = T{ F16_Nguyen_aerodynamics::default_parameters::MAC_m };
        basic_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::wingspan_m] = T{ F16_Nguyen_aerodynamics::default_parameters::wingspan_m };


        // calculate a pitch seed from an aoa seed
        const auto CW = std::cos(tr7::deg2rad(u[input_names::flight_path_angle_deg])) / qS_div_W;

        const auto aoa_point_deg = aoa_wind_point_seed(plant,
                                                       CW,
                                                       q_div_p,
                                                       basic_aeroforce_and_moment_inputs,
                                                       inflags);

        const auto pitch_point_deg = aoa_point_deg + u[input_names::flight_path_angle_deg];


        // calculate a stabilator seed
        auto dhseed_aeroforce_and_moment_inputs = basic_aeroforce_and_moment_inputs;
        dhseed_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::aoa_deg] = aoa_point_deg;
        dhseed_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::dlef_deg] = stationary_lef_schedule(aoa_point_deg, q_div_p, inflags);

        const auto dh_point_deg = stabilator_point_seed(plant, dhseed_aeroforce_and_moment_inputs);


        // solve the trim of two variables (pitch, dh) with the numerical 2d Newton method
        const auto straight_and_level_costs = [&](const auto &pitch_and_dh_deg)
        {
            auto aeroforce_and_moment_inputs = basic_aeroforce_and_moment_inputs;
            aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::aoa_deg] = pitch_and_dh_deg[0] - u[input_names::flight_path_angle_deg];
            aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::dlef_deg] = stationary_lef_schedule(aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::aoa_deg], q_div_p, inflags);
            aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::dh_deg] = pitch_and_dh_deg[1];

            const auto coeffs = plant.aerodynamics().force_and_moment_coefficients(aeroforce_and_moment_inputs);
            const auto &CZtot = coeffs[F16_Nguyen_aerodynamics::force_and_moment_coefficient_names::CZtot];
            const auto &CMtot = coeffs[F16_Nguyen_aerodynamics::force_and_moment_coefficient_names::CMtot];

            return std::array<T, 2>{ qS_div_W * CZtot + std::cos(tr7::deg2rad(pitch_and_dh_deg[0])), CMtot };

        };

        const auto seeds = std::array<T, 2>{ pitch_point_deg, dh_point_deg };
        const auto [pitch_and_dh_deg,
                    costs,
                    exitflag] = tr7::numjac_Newton_method2d(straight_and_level_costs, seeds);

        const auto aero_trim_success = exitflag == tr7::numeric_exitflag::success;
        auto pitch_deg = pitch_and_dh_deg[0];
        auto dh_deg = pitch_and_dh_deg[1];
        auto cost_Fz = costs[0];
        auto cost_My = costs[1];

        if (!aero_trim_success) {
            const auto seed_costs = straight_and_level_costs(seeds);

            if (std::abs(seed_costs[0]) + std::abs(seed_costs[1]) < std::abs(cost_Fz) + std::abs(cost_My)) {
                pitch_deg = seeds[0];
                dh_deg = seeds[1];
                cost_Fz = seed_costs[0];
                cost_My = seed_costs[1];

            }

        }


        // final aero outputs
        const auto aoa_deg = pitch_deg - u[input_names::flight_path_angle_deg];
        const auto dlef_deg = stationary_lef_schedule(aoa_deg, q_div_p, inflags);


        // calculate the necessary thrust for the flight
        auto final_aeroforce_and_moment_inputs = basic_aeroforce_and_moment_inputs;
        final_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::aoa_deg] = aoa_deg;
        final_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::dlef_deg] = dlef_deg;
        final_aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::dh_deg] = dh_deg;

        const auto CXtot = plant.aerodynamics().force_and_moment_coefficients(final_aeroforce_and_moment_inputs)[F16_Nguyen_aerodynamics::force_and_moment_coefficient_names::CXtot];
        const auto necessary_thrust_N = weight_N * (std::sin(tr7::deg2rad(pitch_deg)) - CXtot * qS_div_W);


        // we now call the engine trim
        const auto [engine_trim_success,
                    P3_percent, throttle_percent,
                    cost_Fx] = engine_trim(plant, necessary_thrust_N, TAS_mps / atmos.speed_of_sound_mps, ZP_m);


        // finalize the outputs
        outputs_type<T> trimout{};
        trimout[output_names::cost_Fx] = cost_Fx;
        trimout[output_names::cost_Fz] = cost_Fz;
        trimout[output_names::cost_My] = cost_My;
        trimout[output_names::pitch_deg] = pitch_deg;
        trimout[output_names::dh_deg] = dh_deg;
        trimout[output_names::aoa_deg] = aoa_deg;
        trimout[output_names::dlef_deg] = dlef_deg;
        trimout[output_names::P3_percent] = P3_percent;
        trimout[output_names::throttle_percent] = throttle_percent;

        return std::make_pair(aero_trim_success && engine_trim_success, trimout);

    }


    template<typename F16NguyenPlantType, typename T>
    inline std::tuple<bool, T, T, T> F16_Nguyen_trim::engine_trim(const F16NguyenPlantType &plant,
                                                                  const T &necessary_thrust_N,
                                                                  const T &Mach, const T &ZP_m)
    {
        //
        // x-body engine trim for a steady flight with known
        // thrust. The stationary solution "P3dot = 0" can only
        // be achieved if "P1(throttle) = P2 = P3(thrust)",
        // so we can obtain state "P3" and input "throttle"
        // from the thrust. Returns an std::tuple holding
        // "[trim_success, P3_percent, throttle_percent, trim_cost]".
        //

        // initialize output flag
        auto engine_trim_success = true;


        // create the basic inputs to evaluate the engine dataset
        using engine_inputs_type = F16_Nguyen_engine::inputs_type<T>;
        using engine_input_names = F16_Nguyen_engine::input_names;

        engine_inputs_type basic_engine_inputs{};
        basic_engine_inputs[engine_input_names::ZP_m] = ZP_m;
        basic_engine_inputs[engine_input_names::Mach] = Mach;


        // calculate the engine's regime and "P3 = P3(thrust)"
        const auto engine_dataset_outputs = plant.engine().dataset(basic_engine_inputs);
        const auto &Tidle_N = engine_dataset_outputs[F16_Nguyen_engine::dataset_output_names::Tidle_N];
        const auto &Tmil_N = engine_dataset_outputs[F16_Nguyen_engine::dataset_output_names::Tmil_N];
        const auto &Tmax_N = engine_dataset_outputs[F16_Nguyen_engine::dataset_output_names::Tmax_N];

        const auto P3_high_nondim = (necessary_thrust_N - Tmil_N) / (Tmax_N - Tmil_N) + T{ 1 };
        const auto high_regime = P3_high_nondim >= T{ 1 };

        const auto P3_low_nondim = (necessary_thrust_N - Tidle_N) / (Tmil_N - Tidle_N);
        const auto low_regime = P3_low_nondim <= T{ 1 };

        T P3_percent;
        if ((high_regime && !low_regime) ||
           (high_regime && low_regime && tr7::eq_fp(P3_high_nondim, P3_low_nondim, T{ 1e3 }))) {

            P3_percent = P3_high_nondim * T{ 50 };

        }
        else if (low_regime && !high_regime) {
            P3_percent = P3_low_nondim * T{ 50 };

        }
        else {
            // the engine's P3 regime is UNFEASIBLE, results will be approximate
            engine_trim_success = false;

            const auto in_high = P3_high_nondim - T{ 1 };
            const auto in_low = T{ 1 } - P3_low_nondim;
            P3_percent = in_high > in_low ? P3_high_nondim * T{ 50 } : P3_low_nondim * T{ 50 };

        }

        if (P3_percent > T{ 100 } || P3_percent < T{ 0 }) {
            // the engine's P3 regime is UNFEASIBLE, results will be approximate
            engine_trim_success = false;
            P3_percent = std::clamp(P3_percent, T{ 0 }, T{ 100 });

            const auto &throttle_percent = P3_percent;
            return std::make_tuple(engine_trim_success,
                                   P3_percent, throttle_percent,
                                   std::numeric_limits<T>::infinity());

        }


        // finalize by extracting the throttle from the
        // "P3 = P1(throttle)" identity, solving the curve
        // with "fzero" in the interval "[0, 100]" (i.e.,
        // the domain of the throttle)
        const auto engine_throttle_cost = [&](const auto &throttle_percent)
        {
            auto engine_inputs = basic_engine_inputs;
            engine_inputs[engine_input_names::throttle_percent] = throttle_percent;

            const auto P1_percent = plant.engine().dataset(engine_inputs)[F16_Nguyen_engine::dataset_output_names::P1_percent];

            return P1_percent - P3_percent;

        };

        auto [throttle_percent, cost_engine, exitflag] = tr7::fzero(engine_throttle_cost,
                                                                    T{ 0 }, T{ 100 });

        if (exitflag != tr7::numeric_exitflag::success) {
            engine_trim_success = false;

            const auto seed_cost = engine_throttle_cost(P3_percent);
            if (std::abs(seed_cost) < std::abs(cost_engine)) {
                cost_engine = seed_cost;
                throttle_percent = P3_percent;

            }

        }

        return std::make_tuple(engine_trim_success,
                               P3_percent, throttle_percent,
                               cost_engine);

    }


    template<typename InflagsType, typename T>
    constexpr T F16_Nguyen_trim::stationary_lef_schedule(const T &aoa_deg, const T &q_div_p,
                                                         const InflagsType &inflags)
    {
        //
        // Stationary leading edge flap deflection, from the
        // "Nguyen" reference.
        //

        return static_cast<int>(inflags[inflag_names::disable_lef]) ? T{ 0 } :
                                                                      T{ 1.38 } * aoa_deg - T{ 9.05 } * q_div_p + T{ 1.45 };

    }


    template<typename F16NguyenPlantType,
             typename InflagsType, typename T>
    constexpr T F16_Nguyen_trim::aoa_wind_point_seed(const F16NguyenPlantType &plant,
                                                     const T &CW,
                                                     const T &q_div_p,
                                                     const F16_Nguyen_aerodynamics::force_and_moment_inputs_type<T> &basic_aeroforce_and_moment_inputs,
                                                     const InflagsType &inflags,
                                                     const T &Daoa_numjac_deg,
                                                     std::size_t max_iters_point)
    {
        //
        // Calculates a seed for the trimmed angle-of-attack using
        // the simplest wind-point formula:
        //
        //   CL0 + CLalpha * aoa = CW,
        //
        // where "CW" represents the nondimensional weight,
        //
        //  CW = W * cos(gamma) / cos(mu) / q / S,
        //
        // taking all the signals in "basic_aeroforce_and_moment_inputs"
        // as if they were constant, with the exceptions of the angle of
        // attack and the deflection of the leading edge flap.
        //

        using aeroforce_and_moment_input_names = F16_Nguyen_aerodynamics::force_and_moment_input_names;


        // helper function to calculate CL = CL(aoa)
        const auto CLtot = [&](auto aoa_deg)
        {
            auto aeroforce_and_moment_inputs = basic_aeroforce_and_moment_inputs;
            aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::aoa_deg] = aoa_deg;
            aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::dlef_deg] = stationary_lef_schedule(aoa_deg, q_div_p, inflags);

            const auto coeffs = plant.aerodynamics().force_and_moment_coefficients(aeroforce_and_moment_inputs);

            const auto &CXtot = coeffs[F16_Nguyen_aerodynamics::force_and_moment_coefficient_names::CXtot];
            const auto &CZtot = coeffs[F16_Nguyen_aerodynamics::force_and_moment_coefficient_names::CZtot];
            const auto aoa_rad = tr7::deg2rad(aoa_deg);

            return CXtot * std::sin(aoa_rad) - CZtot * std::cos(aoa_rad);

        };


        // do the thing
        auto aoa_point_deg = T{ 0 };
        for (auto iter_point = 0u; iter_point < max_iters_point; ++iter_point) {
            const auto approx_CLalpha = T{ 0.5 } * (CLtot(aoa_point_deg + Daoa_numjac_deg) - CLtot(aoa_point_deg - Daoa_numjac_deg)) / Daoa_numjac_deg;

            const auto approx_CL0 = CLtot(aoa_point_deg) - approx_CLalpha * aoa_point_deg;

            aoa_point_deg = (CW - approx_CL0) / approx_CLalpha;

        }

        return aoa_point_deg;

    }


    template<typename F16NguyenPlantType, typename T>
    constexpr T F16_Nguyen_trim::stabilator_point_seed(const F16NguyenPlantType &plant,
                                                       const F16_Nguyen_aerodynamics::force_and_moment_inputs_type<T> &basic_aeroforce_and_moment_inputs,
                                                       const T &Ddh_numjac_deg,
                                                       std::size_t max_iters_point)
    {
        //
        // Calculates a seed for the trimmed stabilator using
        // the simplest point formula:
        //
        //   CM0 + CMdh * dh = 0,
        //
        // taking all the signals in "basic_aeroforce_and_moment_inputs"
        // as if they were constant, with the exception of the stabilator
        // deflection.
        //

        using aeroforce_and_moment_input_names = F16_Nguyen_aerodynamics::force_and_moment_input_names;


        // helper function to calculate CM = CM(dh)
        const auto CMtot = [&](auto dh_deg)
        {
            auto aeroforce_and_moment_inputs = basic_aeroforce_and_moment_inputs;
            aeroforce_and_moment_inputs[aeroforce_and_moment_input_names::dh_deg] = dh_deg;

            const auto coeffs = plant.aerodynamics().force_and_moment_coefficients(aeroforce_and_moment_inputs);

            return coeffs[F16_Nguyen_aerodynamics::force_and_moment_coefficient_names::CMtot];

        };


        // do the thing
        auto dh_point_deg = T{ 0 };
        for (auto iter_point = 0u; iter_point < max_iters_point; ++iter_point) {
            const auto approx_CMdh = T{ 0.5 } *(CMtot(dh_point_deg + Ddh_numjac_deg) - CMtot(dh_point_deg - Ddh_numjac_deg)) / Ddh_numjac_deg;

            const auto approx_CM0 = CMtot(dh_point_deg) - approx_CMdh * dh_point_deg;

            dh_point_deg = -approx_CM0 / approx_CMdh;

        }

        return dh_point_deg;

    }


}


#endif
