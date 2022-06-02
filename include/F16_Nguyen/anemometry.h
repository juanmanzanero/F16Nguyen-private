#ifndef F16_NUGYEN_ANEMOMETRY_H
#define F16_NUGYEN_ANEMOMETRY_H
#pragma once


#include "tr7/tr7numeric.h"

#include "F16_Nguyen/isa_atmosphere.h"


//
// Defines useful functions to calculate and convert between
// the typical flight velocities.
//


namespace F16_Nguyen
{
    template<typename T>
    inline T cas(const T &Mach,
                 const T &pressure_Pa,
                 const T &density_SL_kgpm3 = T{ air_density_SL_kgpm3 },
                 const T &pressure_SL_Pa = T{ air_pressure_SL_Pa },
                 const T &gamma = T{ gamma_air })
    {
        //
        // Returns the calibrates airspeed (CAS) given the true Mach number (formed with
        // the true airspeed TAS) and the atmospheric pressure. User can optionally input
        // the sea level air density, pressure and heat capacity ratio. All in SI units.
        //
        // Algorithm:
        //
        //  * The CAS is the speed that the Pitot tube measures. A Pitot tube obtains a 
        //    speed from the dynamic_pressure == impact_pressure = (stagnation_pressure - static_pressure).
        //
        //  * Isentropic flow (with Rankine-Hugoniot for normal shockwaves in supersonic regime) 
        //    formulas yield
        //
        //         impact_pressure = F(Mach, pressure),
        //
        //    and the CAS, by definition, can be obtained solving
        //
        //         impact_pressure = F(CAS / soundspeed_SL, pressure_SL).
        //

        const auto gm1 = gamma - T{ 1 };
        const auto gp1 = gamma + T{ 1 };
        const auto soundspeed2_SL = gamma * pressure_SL_Pa / density_SL_kgpm3;
        const auto Mach2 = Mach * Mach;


        // subsonic impact pressure evaluated at SL and Mach = 1.
        const auto qlimit_SL = pressure_SL_Pa * (tr7::AD_pow(T{ 1 } + T{ 0.5 } * gm1, gamma / gm1) - T{ 1 });


        // evaluate impact pressure
        T impact_pressure{};
        if (Mach <= T{ 1 }) {
            impact_pressure = pressure_Pa * (tr7::AD_pow(T{ 1 } + T{ 0.5 } * gm1 * Mach2, gamma / gm1) - T{ 1 });

        }
        else {
            const auto p_quot_rankine = T{ 1 } + T{ 2 } * gamma / gp1 * (Mach2 - T{ 1 }); // p_postshockwave / p_flight
            const auto p_quot_isentropic = tr7::AD_pow(gp1 * gp1 * Mach2 / (T{ 4 } * gamma * Mach2 - T{ 2 } * gm1), gamma / gm1); // p_stagnation / p_postshockwave

            impact_pressure = pressure_Pa * (p_quot_isentropic * p_quot_rankine - T{ 1 });

        }


        // evaluate CAS
        if (impact_pressure <= qlimit_SL) {
            return tr7::AD_sqrt(T{ 2 } * soundspeed2_SL / gm1 * (tr7::AD_pow(impact_pressure / pressure_SL_Pa + T{ 1 }, gm1 / gamma) - T{ 1 }));

        }
        else {
            // the following equation must be solved:
            // k * d = y / (1 - b / y)^n
            // y = calibrated_Mach^2 = (cas_mps / soundspeed_SL)^2
            // d = impact_pressure / pressure_SL + 1
            // k, b, n = f(gamma)

            const auto d = impact_pressure / pressure_SL_Pa + T{ 1 };
            const auto n = T{ 1 } / gm1;
            const auto k = T{ 1 } / tr7::AD_pow(T{ 0.5 } * gp1, gamma / gm1) / tr7::AD_pow(T{ 0.5 } * gp1 / gamma, n); // 1. / (1.2^3.5 * (6 / 7)^2.5)
            const auto b = T{ 0.5 } * gm1 / gamma;
            const auto kd = k * d;


            // solve the equation numerically
#ifdef TR7_WITH_AD
            const auto n_numerical = tr7::AD_force_scalar_value(n);
            const auto b_numerical = tr7::AD_force_scalar_value(b);
            const auto kd_numerical = tr7::AD_force_scalar_value(kd);

#else
            const auto &n_numerical = n;
            const auto &b_numerical = b;
            const auto &kd_numerical = kd;
#endif
            const auto kdnb_numerical = kd_numerical * n_numerical * b_numerical;

            constexpr auto max_evals = std::size_t{ 30u };
            constexpr auto tol = decltype(n_numerical){ 1e-8 };
            constexpr auto calibrated_Mach_seed = decltype(n_numerical){ 1 };
            constexpr auto seed = calibrated_Mach_seed * calibrated_Mach_seed;

            const auto F = [=](const auto &y_numerical)
            {
                return kd_numerical * std::pow(decltype(y_numerical){ 1 } - b_numerical / y_numerical,
                                               n_numerical) - y_numerical;

            };


            // Newton method
            const auto dF_dy = [=](const auto &y_numerical)
            {
                return kdnb_numerical * std::pow(decltype(y_numerical){ 1 } - b_numerical / y_numerical,
                                                 n_numerical - decltype(y_numerical){ 1 }) / y_numerical / y_numerical -
                                                 decltype(y_numerical){ 1 };

            };

            const auto [y_numerical, __ignore, exitflag] = tr7::Newton_method(F, dF_dy, seed, decltype(n_numerical){ 1 },
                                                                              max_evals, tol, tol); (void)__ignore;


            T y{ y_numerical };
            if constexpr (tr7::is_AD_scalar<T>::value) {
                if (exitflag == tr7::numeric_exitflag::success) {
                    // apply the derivatives of "y" with respect to "kd", "b" and "n"
                    const auto f = decltype(y_numerical){ 1 } - b_numerical / y_numerical;

                    const auto dF_dkd = std::pow(f, n_numerical);

                    const auto dF_db = -kd_numerical * n_numerical * dF_dkd / (f * y_numerical);

                    const auto dF_dn = kd_numerical * dF_dkd * tr7::AD_log(f);

                    y -= (dF_dn * (n - n_numerical) +
                          dF_dkd * (kd - kd_numerical) +
                          dF_db * (b - b_numerical)) / dF_dy(y_numerical);

                }

            }

            if (exitflag != tr7::numeric_exitflag::success) {
                // revert to a regression (I found it in the EFA code...)
                y = T{ 0.41726 } + T{ 0.7767 } * (d - T{ 1 }) - T{ 0.0989 } / (d - T{ 1 });

            }

            return tr7::AD_sqrt(soundspeed2_SL * y);

        }

    }


    template<typename T>
    constexpr pure_function T subsonic_cas(const T &Mach,
                                           const T &pressure_Pa,
                                           const T &density_SL_kgpm3 = T{ air_density_SL_kgpm3 },
                                           const T &pressure_SL_Pa = T{ air_pressure_SL_Pa },
                                           const T &gamma = T{ gamma_air })
    {
        //
        // Returns the calibrated airspeed (CAS) for subsonic Mach
        // numbers (formed with the true airspeed TAS), atmospheric
        // pressure, and optionally, the sea-level air density,
        // pressure, and specific heat ratio. All in SI units.
        //
        // NOTE: the function ALWAYS applies the subsonic expressions for
        // the CAS, irrespectively of the input values, like for example
        // the boom simulator Simulink model does.
        //

        // helper ideal gas constants
        constexpr_variable auto gm1_div_gamma = (gamma - T{ 1 }) / gamma;

        // subsonic impact pressure & CAS formulas
        const auto impact_pressure = pressure_Pa *
                                     (tr7::AD_pow(T{ 1 } +
                                                  T{ 0.5 } * gamma * gm1_div_gamma * Mach * Mach,
                                                  T{ 1 } / gm1_div_gamma) -
                                      T{ 1 });

        return tr7::AD_sqrt(T{ 2 } * pressure_SL_Pa / density_SL_kgpm3 / gm1_div_gamma *
                            (tr7::AD_pow(impact_pressure / pressure_SL_Pa + T{ 1 }, gm1_div_gamma) - T{ 1 }));

    }


    template<typename T>
    constexpr T density_sigma(const T &geopotential_altitude_m,
                              const T &temperature_offset_K = T{ 0 })
    {
        //
        // Returns the quotient (density / density_SL), given the
        // geopotential altitude and an optional temperature
        // offset. All in SI units.
        //

        static_assert(!tr7::is_AD_scalar<T>::value, "F16_Nguyen::density_sigma: not prepared for AD...");


        // zero altitude solution
        if (tr7::eq_fp(geopotential_altitude_m, T{ 0 })) {
            return T{ 1 };
        }


        // do the thing
        const auto atmos = isa_atmosphere(geopotential_altitude_m, temperature_offset_K);
        if (tr7::eq_fp(temperature_offset_K, T{ 0 })) {
            return atmos.air_density_kgpm3 / T{ air_density_SL_kgpm3 };

        }
        else {
            const auto atmos_SL = isa_atmosphere(T{ 0 }, temperature_offset_K);
            return atmos.air_density_kgpm3 / atmos_SL.air_density_kgpm3;

        }

    }


    template<typename T>
    constexpr T tas2eas(const T &TAS_mps, const T &geopotential_altitude_m, const T &temperature_offset_K = T{ 0 })
    {
        //
        // EAS (equivalent airspeed) from TAS (true airspeed)
        // and geopotential altitude, with optional temperature
        // offset. All in SI units.
        //

        static_assert(!tr7::is_AD_scalar<T>::value, "F16_Nguyen::tas2eas: not prepared for AD...");


        // zero altitude solution
        if (tr7::eq_fp(geopotential_altitude_m, T{ 0 })) {
            return TAS_mps;
        }


        // do the thing
        return TAS_mps * std::sqrt(density_sigma(geopotential_altitude_m, temperature_offset_K));

    }


    template<typename T>
    constexpr T eas2tas(const T &EAS_mps, const T &geopotential_altitude_m, const T &temperature_offset_K = T{ 0 })
    {
        //
        // TAS (true airspeed) from EAS (equivalent airspeed)
        // and geopotential altitude, with optional temperature
        // offset. All in SI units.
        //

        static_assert(!tr7::is_AD_scalar<T>::value, "F16_Nguyen::eas2tas: not prepared for AD...");


        // zero altitude solution
        if (tr7::eq_fp(geopotential_altitude_m, T{ 0 })) {
            return EAS_mps;
        }


        // do the thing
        return EAS_mps / std::sqrt(density_sigma(geopotential_altitude_m, temperature_offset_K));

    }


    template<typename T>
    inline T tas2cas(const T &TAS_mps, const T &geopotential_altitude_m, const T &temperature_offset_K = T{ 0 })
    {
        //
        // CAS (calibrated airspeed) from TAS (true airspeed) and geopotential altitude, 
        // with optional temperature offset. All in SI units.
        //

        static_assert(!tr7::is_AD_scalar<T>::value, "F16_Nguyen::tas2cas: not prepared for AD...");


        // zero altitude solution
        if (tr7::eq_fp(geopotential_altitude_m, T{ 0 })) {
            return TAS_mps;
        }


        // do the thing
        const auto atmos = isa_atmosphere(geopotential_altitude_m, temperature_offset_K);
        if (tr7::eq_fp(temperature_offset_K, T{ 0 })) {
            return cas(TAS_mps / atmos.speed_of_sound_mps, atmos.air_pressure_Pa);

        }
        else {
            const auto atmos_SL = isa_atmosphere(T{ 0 }, temperature_offset_K);
            return cas(TAS_mps / atmos.speed_of_sound_mps, atmos.air_pressure_Pa,
                       atmos_SL.air_density_kgpm3, atmos_SL.air_pressure_Pa);

        }

    }


    template<typename T>
    inline T cas2tas(const T &CAS_mps, const T &geopotential_altitude_m, const T &temperature_offset_K = T{ 0 })
    {
        //
        // TAS (true airspeed) from CAS (calibrated airspeed)
        // and geopotential altitude, with optional temperature
        // offset. All in SI units.
        //

        static_assert(!tr7::is_AD_scalar<T>::value, "F16_Nguyen::cas2tas: not prepared for AD...");


        // zero altitude solution
        if (tr7::eq_fp(geopotential_altitude_m, T{ 0 })) {
            return CAS_mps;
        }


        // do the thing
        constexpr auto max_evals = std::size_t{ 30u };
        constexpr auto tol = 1e-8;
        const auto seed = eas2tas(CAS_mps, geopotential_altitude_m, temperature_offset_K);

        const auto F = [=](const auto &TAS_mps) { return tas2cas(TAS_mps, geopotential_altitude_m, temperature_offset_K) - CAS_mps; };

        const auto [TAS_mps, __ignore, exitflag] = tr7::fzero(F, seed,
                                                              max_evals,
                                                              tol, tol); (void)__ignore;

        if (exitflag != tr7::numeric_exitflag::success) {
            // revert to the seed (equal to TAS if we took the input as being EAS instead of CAS)
            return seed;

        }
        else {
            return TAS_mps;

        }

    }


    template<typename T>
    inline T eas2cas(const T &EAS_mps, const T &geopotential_altitude_m, const T &temperature_offset_K = T{ 0 })
    {
        //
        // CAS (calibrated airspeed) from EAS (equivalent airspeed)
        // and geopotential altitude, with optional temperature
        // offset. All in SI units.
        //

        static_assert(!tr7::is_AD_scalar<T>::value, "F16_Nguyen::eas2cas: not prepared for AD...");


        // zero altitude solution
        if (tr7::eq_fp(geopotential_altitude_m, T{ 0 })) {
            return EAS_mps;
        }


        // do the thing
        const auto atmos = isa_atmosphere(geopotential_altitude_m, temperature_offset_K);
        if (tr7::eq_fp(temperature_offset_K, T{ 0 })) {
            return cas(EAS_mps / std::sqrt(atmos.air_density_kgpm3 / T{ air_density_SL_kgpm3 }) / atmos.speed_of_sound_mps,
                       atmos.air_pressure_Pa);

        }
        else {
            const auto atmos_SL = isa_atmosphere(T{ 0 }, temperature_offset_K);
            return cas(EAS_mps / std::sqrt(atmos.air_density_kgpm3 / atmos_SL.air_density_kgpm3) / atmos.speed_of_sound_mps,
                       atmos.air_pressure_Pa,
                       atmos_SL.air_density_kgpm3, atmos_SL.air_pressure_Pa);

        }

    }


    template<typename T>
    inline T cas2eas(const T &CAS_mps, const T &geopotential_altitude_m, const T &temperature_offset_K = T{ 0 })
    {
        //
        // EAS (equivalent airspeed) from CAS (calibrated airspeed)
        // and geopotential altitude, with optional temperature offset.
        // All in SI units.
        //

        static_assert(!tr7::is_AD_scalar<T>::value, "F16_Nguyen::cas2eas: not prepared for AD...");

        // zero altitude solution
        if (tr7::eq_fp(geopotential_altitude_m, T{ 0 })) {
            return CAS_mps;
        }


        // do the thing
        constexpr auto max_evals = std::size_t{ 30u };
        constexpr auto tol = 1e-8;
        const auto seed = CAS_mps;

        const auto F = [=](const auto &EAS_mps) { return eas2cas(EAS_mps, geopotential_altitude_m, temperature_offset_K) - CAS_mps; };

        const auto [EAS_mps, __ignore, exitflag] = tr7::fzero(F, seed,
                                                              max_evals,
                                                              tol, tol); (void)__ignore;

        if (exitflag != tr7::numeric_exitflag::success) {
            // revert to the seed (the input CAS)
            return seed;

        }
        else {
            return EAS_mps;

        }

    }


    template<typename T>
    constexpr pure_function T subsonic_cas2Mach(const T &CAS_mps,
                                                const T &pressure_Pa,
                                                const T &temperature_SL_K = T{ air_temperature_SL_K },
                                                const T &pressure_SL_Pa = T{ air_pressure_SL_Pa },
                                                const T &gamma = T{ gamma_air },
                                                const T &Rg = T{ Rg_air_JpkgpK })
    {
        //
        // Returns the Mach number (formed with the true airspeed) from
        // input CAS (calibrated airspeed), pressure, and optional sea-level
        // temperature, pressure, specific heat ratio and ideal gas constant.
        // All in SI units.
        //
        // NOTE: the function ALWAYS applies the subsonic expressions for
        // the CAS, irrespectively of the input values, like for example
        // the boom simulator Simulink model does.
        //

        // helper ideal gas constants
        constexpr_variable auto gm1_div_gamma = (gamma - T{ 1 }) / gamma;

        // subsonic Mach formula
        return tr7::AD_sqrt(T{ 2 } / gamma / gm1_div_gamma * (tr7::AD_pow(T{ 1 } +
                                                                          pressure_SL_Pa / pressure_Pa * (tr7::AD_pow(T{ 1 } +
                                                                                                                      T{ 0.5 } * gm1_div_gamma * CAS_mps * CAS_mps / temperature_SL_K / Rg,
                                                                                                                      T{ 1 } / gm1_div_gamma) - T{ 1 }),
                                                                          gm1_div_gamma) - T{ 1 }));

    }


}


#endif