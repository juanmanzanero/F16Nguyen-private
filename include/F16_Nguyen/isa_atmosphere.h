#ifndef F16_NUGYEN_ISA_ATMOSPHERE_H
#define F16_NUGYEN_ISA_ATMOSPHERE_H
#pragma once


#include "tr7/tr7math.h"


//
// Defines constants and functions that
// deal with the ISA standard atmosphere
// model.
//


namespace F16_Nguyen
{
    //
    // Some common constants of the atmosphere model, we'll
    // declare them as doubles.
    //

    inline constexpr double g0_mps2 = 9.80665; // standard acceleration of gravity [m / s^2]
    inline constexpr double Rg_air_JpkgpK = 287.0531; // atmospheric air ideal gas constant [J / kg / K]
    inline constexpr double gamma_air = 1.4; // atmospheric air ideal gas specific heat ratio
    inline constexpr double air_pressure_SL_Pa = 101325; // standard sea-level atmospheric pressure [Pa]
    inline constexpr double air_temperature_SL_K = 288.15; // standard sea-level atmospheric temperature [K]
    inline constexpr double air_density_SL_kgpm3 = 1.224999036598556; // standard sea-level atmospheric density [kg / m^3] (p / R * T)


    //
    // The aggregate of results that the atmosphere
    // functions return.
    //

    template<typename T>
    struct atmosphere_results
    {
        T air_density_kgpm3;
        T speed_of_sound_mps;
        T air_kinematic_viscosity_m2ps;
        T air_temperature_K;
        T air_pressure_Pa;

    };


    //
    // Atmosphere model functions.
    //

    template<typename T>
    constexpr pure_function T geopotential_altitude(const T &geometric_altitude_m)
    {
        //
        // Returns the Geopotential altitude [m] calculated
        // from the input geometric altitude [m].
        //

        constexpr_variable auto Earth_equatorial_radius_m = T{ 6378137 }; // m
        constexpr_variable auto Earth_polar_radius_m = T{ 63567523 }; // m
        constexpr_variable auto Earth_mean_radius_m = (T{ 2 } * Earth_equatorial_radius_m + Earth_polar_radius_m) / T{ 3 }; // m
 
        return geometric_altitude_m * Earth_mean_radius_m / (geometric_altitude_m + Earth_mean_radius_m);

    }


    template<typename T>
    constexpr pure_function T Sutherland_law(const T &temperature_K)
    {
        //
        // Sutherland's viscosity law for the dynamic viscosity
        // of a gas ("mu"), in SI units [kg / m / s].
        //

        constexpr_variable auto S = T{ 110.4 }; // Sutherland's temperature [K]

        //constexpr_variable auto mu0 = T{ 1.716e-5 }; // kg / m / s
        //constexpr_variable auto T_S = T{ 273.15 }; // K
        constexpr_variable auto C1 = T{ 1.457932654517625e-06 }; // kg / m / s / K^(1 / 2), = mu0 * (T_S + S) / T_S / sqrt(T_S);

        return C1 * temperature_K * tr7::AD_sqrt(temperature_K) / (temperature_K + S);

    }


    template<typename T>
    constexpr atmosphere_results<T> isa_atmosphere(const T &geopotential_altitude_m,
                                                   const T &temperature_offset_K = T{ 0 })
    {
        //
        // Calculates results for the standard ISA model.
        //

        // helper model constants
        constexpr_variable auto g0_div_R = T{ g0_mps2 / Rg_air_JpkgpK };

        // lapse rates [K / m]
        constexpr auto num_atmospheric_layers = std::size_t{ 8u };
        constexpr_variable auto K = std::array<T, num_atmospheric_layers>{
                                        T{ -0.0065 }, // Troposphere
                                        T{ 0 },       // Tropopause
                                        T{ 0.001 },   // Stratosphere1
                                        T{ 0.0028 },  // Stratosphere2
                                        T{ 0 },       // Stratopause
                                        T{ -0.0028 }, // Mesosphere1
                                        T{ -0.002 },  // Mesosphere2
                                        T{ 0. } };    // Mesopause

        // base temperatures [K]
        constexpr_variable auto TH = std::array<T, num_atmospheric_layers>{
                                         T{ air_temperature_SL_K }, // Troposphere
                                         T{ 216.65 },               // Tropopause
                                         T{ 216.65 },               // Stratosphere1
                                         T{ 228.65 },               // Stratosphere2
                                         T{ 270.65 },               // Stratopause
                                         T{ 270.65 },               // Mesosphere1
                                         T{ 214.65 },               // Mesosphere2
                                         T{ 186.94590831019 } };    // Mesopause

        // base geopotential altitudes [m]
        constexpr_variable auto H = std::array<T, num_atmospheric_layers>{
                                       T{ 0 } ,                 // Troposphere
                                       T{ 11000 },              // Tropopause
                                       T{ 20000 },              // Stratosphere1
                                       T{ 32000 },              // Stratosphere2
                                       T{ 47000 },              // Stratopause
                                       T{ 51000 },              // Mesosphere1
                                       T{ 71000 },              // Mesosphere2
                                       T{ 84852.0458449057 } }; // Mesopause

        // base pressures [Pa]
        constexpr_variable auto P = std::array<T, num_atmospheric_layers>{
                                        T{ air_pressure_SL_Pa },  // Troposphere
                                        T{ 22632.0400950078 },    // Tropopause
                                        T{ 5474.87742428105 },    // Stratosphere1
                                        T{ 868.015776620216 },    // Stratosphere2
                                        T{ 110.90577336731 },     // Stratopause
                                        T{ 66.9385281211797 },    // Mesosphere1
                                        T{ 3.9563921603966 },     // Mesosphere2
                                        T{ 0.373377173762337 } }; // Mesopause


        // find the layer we're in
        const auto layer = std::distance(H.cbegin(),
                                         std::lower_bound(std::next(H.cbegin()), H.cend(),
                                         geopotential_altitude_m)) - 1u;


        // calculate the air's temperature & pressure
        T temperature_K;
        T pressure_Pa;
        if (!tr7::eq_fp(K[layer], T{ 0 }, T{ 1e3 })) {
            // linear temperature layer
            const auto theta = T{ 1 } + K[layer] * (geopotential_altitude_m - H[layer]) / TH[layer];
            temperature_K = theta * TH[layer] + temperature_offset_K;
            pressure_Pa = P[layer] * tr7::AD_pow(theta, -g0_div_R / K[layer]);

        }
        else {
            // constant temperature layer
            temperature_K = TH[layer] + temperature_offset_K;
            pressure_Pa = P[layer] * tr7::AD_exp(-g0_div_R * (geopotential_altitude_m - H[layer]) / TH[layer]);

        }

        const auto RT = T{ Rg_air_JpkgpK } * temperature_K;
        const auto density_kgpm3 = pressure_Pa / RT;


        // gather outputs
        return atmosphere_results<T>{ density_kgpm3,
                                      tr7::AD_sqrt(T{ gamma_air } * RT),
                                      Sutherland_law(temperature_K) / density_kgpm3,
                                      temperature_K,
                                      pressure_Pa };

    }


    template<typename T>
    constexpr atmosphere_results<T> simplified_atmosphere_model(const T &geopotential_altitude_m,
                                                                const T &temperature_offset_K = T{ 0 })
    {
        //
        // Returns results for the atmospheric air, mimicking the
        // "Atmosphere model" Simulink block. The Simulink block's mask
        // claims that it implements the ISA standard model for altitudes
        // from 0 to 20 km, but it does not have the tropopause, and does
        // a strange mix with the stratosphere...
        //

        // helper model constants
        constexpr_variable auto p0_troposphere = T{ air_pressure_SL_Pa }; // base pressure (standard sea-level atmospheric pressure) [Pa]
        constexpr_variable auto T0_troposphere = T{ air_temperature_SL_K }; // base temp. (standard sea-level atmospheric temperature) [K]
        constexpr_variable auto rho0_troposphere = T{ air_density_SL_kgpm3 }; // base density [kg / m^3]
        constexpr_variable auto L0_troposphere = T{ 0.0065 }; // lapse rate [K / m]
        constexpr_variable auto H0_troposphere = T{ 0 }; // troposphere base geop. alt. [m]
        constexpr_variable auto H0_tropopause = T{ 11000 }; // tropopause base geop. alt. [m]
        constexpr_variable auto H0_stratosphere = T{ 20000 }; // stratosphere base geop. alt. [m]


        // limit the input geopotential altitude to the troposphere
        const auto tropoaltitude_m = tr7::AD_clamp(geopotential_altitude_m,
                                                   H0_troposphere, H0_tropopause);


        // temperature
        const auto T_troposphere = T0_troposphere - L0_troposphere * tropoaltitude_m;
        const auto temperature_K = T_troposphere + temperature_offset_K;
        const auto RT = T{ Rg_air_JpkgpK } * temperature_K;


        // pressure & density in the troposphere
        const auto tropomodel = tr7::AD_pow(T_troposphere / T0_troposphere,
                                            T{ g0_mps2 } / (L0_troposphere * T{ Rg_air_JpkgpK }));

        const auto pressure_troposphere_Pa = p0_troposphere * tropomodel;
        const auto density_troposphere_kgpm3 = pressure_troposphere_Pa / RT;


        // the Simulink block now calculates the stratosphere values and a final
        // pressure & density (wtf ??)
        const auto stratoaltitude_m = tr7::AD_clamp(H0_tropopause - geopotential_altitude_m,
                                                    H0_tropopause - H0_stratosphere, H0_troposphere);
        const auto stratomodel = tr7::AD_exp(stratoaltitude_m * T{ g0_mps2 } /
                                             (T{ Rg_air_JpkgpK } * T_troposphere));

        const auto pressure_Pa = stratomodel * pressure_troposphere_Pa;
        const auto density_kgpm3 = stratomodel * density_troposphere_kgpm3;


        // gather outputs
        return atmosphere_results<T>{ density_kgpm3,
                                      tr7::AD_sqrt(T{ gamma_air } * RT),
                                      Sutherland_law(temperature_K) / density_kgpm3,
                                      temperature_K,
                                      pressure_Pa };

    }


}


#endif