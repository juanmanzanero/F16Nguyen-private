#ifndef F16_NGUYEN_F16_NGUYEN_AIR_KINEMATICS_H
#define F16_NGUYEN_F16_NGUYEN_AIR_KINEMATICS_H
#pragma once


#include "F16_Nguyen/modelglue.h"
#include "F16_Nguyen/algebra3d.h"
#include "F16_Nguyen/isa_atmosphere.h"
#include "F16_Nguyen/anemometry.h"
#include "F16_Nguyen/F16_Nguyen_flat_Earth_body_dynamics.h"


//
// Defines the "F16_Nguyen_air_kinematics" class,
// that calculates aero-related variables to use
// as inputs for the evaluation of aerodynamic forces
// & moments.
//


namespace F16_Nguyen
{
    class F16_Nguyen_air_kinematics
    {
        //
        // Calculates aero-related variables to
        // be used as inputs to the aerodynamic
        // forces & moments.
        //

    public:

        DECLARE_MODEL_OUTPUTS(dynamic_pressure_Pa,
                              Mach, KCAS, KEAS, ZP_ft,
                              VTASx_mps, VTASy_mps, VTASz_mps,
                              TAS_mps, aoa_deg, aos_deg,
                              true_heading_angle_deg, flight_path_angle_deg, flight_roll_angle_deg,
                              pw_radps, qw_radps, rw_radps, // body angular velocity, but projected onto the wind-axes frame
                              paero_radps, qaero_radps, raero_radps, // body "aerodynamic" angular velocity, projected onto the body-axes frame (it's the one we should take as input to calculate the aerodynamic forces & moments)
                              air_density_kgpm3, sounspeed_mps,
                              air_kinematic_viscosity_m2ps, air_temperature_K, air_pressure_Pa)


        template<typename T>
        static constexpr outputs_type<T> outputs(const T &altitude_m,
                                                 const T &VGSx_mps, const T &VGSy_mps, const T &VGSz_mps,
                                                 const T &p_radps, const T &q_radps, const T &r_radps,
                                                 const T &qw_body2Earth, const T &qx_body2Earth, const T &qy_body2Earth, const T &qz_body2Earth,
                                                 const T &DISA_K,
                                                 const T &VWx_bodyaxes_m = T{ 0 },
                                                 const T &VWy_bodyaxes_m = T{ 0 },
                                                 const T &VWz_bodyaxes_m = T{ 0 },
                                                 const T &pgust_bodyaxes_radps = T{ 0 },
                                                 const T &qgust_bodyaxes_radps = T{ 0 },
                                                 const T &rgust_bodyaxes_radps = T{ 0 });

    };


    //
    // "F16_Nguyen_air_kinematics" impl.
    //

    template<typename T>
    constexpr auto F16_Nguyen_air_kinematics::outputs(const T &altitude_m,
                                                      const T &VGSx_mps, const T &VGSy_mps, const T &VGSz_mps,
                                                      const T &p_radps, const T &q_radps, const T &r_radps,
                                                      const T &qw_body2Earth, const T &qx_body2Earth, const T &qy_body2Earth, const T &qz_body2Earth,
                                                      const T &DISA_K,
                                                      const T &VWx_bodyaxes_m,
                                                      const T &VWy_bodyaxes_m,
                                                      const T &VWz_bodyaxes_m,
                                                      const T &pgust_bodyaxes_radps,
                                                      const T &qgust_bodyaxes_radps,
                                                      const T &rgust_bodyaxes_radps) -> outputs_type<T>
    {
        //
        // Evaluates the model's output equation, i.e.,
        // "outputs = G(state, inputs)", and returns the
        // model's "y".
        //

        // evaluate the ISA atmosphere model
        constexpr auto geopotential_altitude_formula = false;
        T ZP_m;
        if constexpr (geopotential_altitude_formula) {
            ZP_m = geopotential_altitude(altitude_m);
        }
        else {
            ZP_m = altitude_m;
        }

        const auto atmos = isa_atmosphere(ZP_m, DISA_K);


        // calculate the true airspeed vector
        const auto VTASx_mps = VGSx_mps - VWx_bodyaxes_m;
        const auto VTASy_mps = VGSy_mps - VWy_bodyaxes_m;
        const auto VTASz_mps = VGSz_mps - VWz_bodyaxes_m;


        // calculate the anemometry variables
        const auto TASsqr = VTASx_mps * VTASx_mps +
                            VTASy_mps * VTASy_mps +
                            VTASz_mps * VTASz_mps;

        const auto TAS_mps = tr7::AD_sqrt(TASsqr);
        const auto Mach = TAS_mps / atmos.speed_of_sound_mps;
        const auto dynamic_pressure_Pa = T{ 0.5 } * atmos.air_density_kgpm3 * TASsqr;

        T KEAS{};
        T KCAS{};
        if constexpr (tr7::is_AD_scalar<T>::value) {
            const auto atmos_SL = isa_atmosphere(T{ 0 }, DISA_K);
            KEAS = tr7::mps2kn(TAS_mps * tr7::AD_sqrt(atmos.air_density_kgpm3 / atmos_SL.air_density_kgpm3));
            KCAS = tr7::mps2kn(cas(Mach, atmos.air_pressure_Pa, atmos_SL.air_density_kgpm3, atmos_SL.air_pressure_Pa));

        }
        else {
            if (tr7::eq_fp(ZP_m, T{ 0 })) {
                KEAS = tr7::mps2kn(TAS_mps);
                KCAS = KEAS;

            }
            else if (tr7::eq_fp(DISA_K, T{ 0 })) {
                KEAS = tr7::mps2kn(TAS_mps * tr7::AD_sqrt(atmos.air_density_kgpm3 / T{ air_density_SL_kgpm3 }));
                KCAS = tr7::mps2kn(cas(Mach, atmos.air_pressure_Pa));

            }
            else {
                const auto atmos_SL = isa_atmosphere(T{ 0 }, DISA_K);
                KEAS = tr7::mps2kn(TAS_mps * tr7::AD_sqrt(atmos.air_density_kgpm3 / atmos_SL.air_density_kgpm3));
                KCAS = tr7::mps2kn(cas(Mach, atmos.air_pressure_Pa, atmos_SL.air_density_kgpm3, atmos_SL.air_pressure_Pa));

            }

        }


        // calculate the angles of attack and
        // sideslip, and he flight path angles
        const auto [aoa_rad,
                    aos_rad] = velocity_angles(VTASx_mps, VTASy_mps, VTASz_mps);

        const auto T_body2Earth = quat2rotmat(qw_body2Earth, qx_body2Earth, qy_body2Earth, qz_body2Earth);

        const auto T_body2wind = velocity_angles2rotmat(aoa_rad, aos_rad);

        const auto T_wind2Earth = matrix3x3_times_matrix3x3(T_body2Earth,
                                                            matrix3x3_transpose(T_body2wind));

        const auto [true_heading_angle_rad,
                    flight_path_angle_rad,
                    flight_roll_angle_rad] = rotmat2ea(T_wind2Earth);


        // project the body angular velocities onto the wind-axes frame
        const auto [pw_radps,
                    qw_radps,
                    rw_radps] = matrix3x3_times_vector3d(T_body2wind,
                                                         p_radps, q_radps, r_radps);


        // calculate the "aerodynamic" body-angular velocity
        const auto paero_radps = p_radps - pgust_bodyaxes_radps;
        const auto qaero_radps = q_radps - qgust_bodyaxes_radps;
        const auto raero_radps = r_radps - rgust_bodyaxes_radps;


        // return the model's outputs
        return outputs_type<T>{ dynamic_pressure_Pa,
                                Mach, KCAS, KEAS, tr7::m2ft(ZP_m),
                                VTASx_mps, VTASy_mps, VTASz_mps,
                                TAS_mps, tr7::rad2deg(aoa_rad), tr7::rad2deg(aos_rad),
                                tr7::rad2deg(true_heading_angle_rad), tr7::rad2deg(flight_path_angle_rad), tr7::rad2deg(flight_roll_angle_rad),
                                pw_radps, qw_radps, rw_radps,
                                paero_radps, qaero_radps, raero_radps,
                                atmos.air_density_kgpm3, atmos.speed_of_sound_mps,
                                atmos.air_kinematic_viscosity_m2ps, atmos.air_temperature_K, atmos.air_pressure_Pa };

   }


}


#endif