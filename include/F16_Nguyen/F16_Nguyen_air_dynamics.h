#ifndef F16_NGUYEN_F16_NGUYEN_AIR_DYNAMICS_H
#define F16_NGUYEN_F16_NGUYEN_AIR_DYNAMICS_H
#pragma once


#include "F16_Nguyen/modelglue.h"
#include "F16_Nguyen/algebra3d.h"


//
// Defines the "F16_Nguyen_air_dynamics" class,
// that calculates aero-related variables that
// take the acceleration into account (this is
// what distinguishes it from the model in
// "F16_Nguyen_air_kinematics.h").
//


namespace F16_Nguyen
{
    class F16_Nguyen_air_dynamics
    {
        //
        // A class to calculate some useful
        // aerodynamic outputs after we know
        // the body's accelerations.
        //

    public:

        DECLARE_MODEL_OUTPUTS(dVTASdtx_mps2, dVTASdty_mps2, dVTASdtz_mps2, // total time derivative components of the true airspeed vector
                              VTASxdot_mps2, VTASydot_mps2, VTASzdot_mps2, // pseudo-accelerations of the true airspeed vector
                              TASdot_mps2, aoadot_degps, aosdot_degps,
                              Nwx_g, Nwy_g, Nwz_g,  // body load factors, but projected onto the wind-axes frame
                              p_windREarth_radps, q_windREarth_radps, r_windREarth_radps, // angular velocities of the wind axes frame relative to the Earth axes one, projected onto the former
                              true_heading_angledot_degps, flight_path_angledot_degps, flight_roll_angledot_degps,
                              air_turn_radius_m)


        template<typename T>
        static constexpr outputs_type<T> outputs(const T &VGSdotx_mps, const T &VGSdoty_mps, const T &VGSdotz_mps,
                                                 const T &VTASx_mps, const T &VTASy_mps, const T &VTASz_mps,
                                                 const T &TAS_mps, const T &aoa_deg, const T &aos_deg,
                                                 const T &true_heading_angle_deg, const T &flight_path_angle_deg, const T &flight_roll_angle_deg,
                                                 const T &Nx_g, const T &Ny_g, const T &Nz_g,
                                                 const T &p_radps, const T &q_radps, const T &r_radps,
                                                 const T &VWdotx_bodyaxes_m = T{ 0 },
                                                 const T &VWdoty_bodyaxes_m = T{ 0 },
                                                 const T &VWdotz_bodyaxes_m = T{ 0 });

    };


    //
    // "F16_Nguyen_air_dynamics" impl.
    //

    template<typename T>
    constexpr auto F16_Nguyen_air_dynamics::outputs(const T &VGSdotx_mps, const T &VGSdoty_mps, const T &VGSdotz_mps,
                                                    const T &VTASx_mps, const T &VTASy_mps, const T &VTASz_mps,
                                                    const T &TAS_mps, const T &aoa_deg, const T &aos_deg,
                                                    const T &true_heading_angle_deg, const T &flight_path_angle_deg, const T &flight_roll_angle_deg,
                                                    const T &Nx_g, const T &Ny_g, const T &Nz_g,
                                                    const T &p_radps, const T &q_radps, const T &r_radps,
                                                    const T &VWdotx_bodyaxes_m,
                                                    const T &VWdoty_bodyaxes_m,
                                                    const T &VWdotz_bodyaxes_m) -> outputs_type<T>
    {
        //
        // Evaluates the model's output equation, i.e.,
        // "outputs = G(state, inputs)", and returns the
        // model's "y".
        //

        // calculate the body-to-wind rotation matrix
        const auto aoa_rad = tr7::deg2rad(aoa_deg);
        const auto aos_rad = tr7::deg2rad(aos_deg);
        const auto T_body2wind = velocity_angles2rotmat(aoa_rad,
                                                        aos_rad);


        // calculate the pseudo-accelerations of the
        // true airspeed vector
        const auto VTASxdot_mps2 = VGSdotx_mps - VWdotx_bodyaxes_m;
        const auto VTASydot_mps2 = VGSdoty_mps - VWdoty_bodyaxes_m;
        const auto VTASzdot_mps2 = VGSdotz_mps - VWdotz_bodyaxes_m;


        // calculate the total time derivative of the
        // true airspeed vector
        const auto dVTASdt_mps2 = pseudoacceleration2acceleration(VTASxdot_mps2, VTASydot_mps2, VTASzdot_mps2,
                                                                  VTASx_mps, VTASy_mps, VTASz_mps,
                                                                  p_radps, q_radps, r_radps);


        // calculate the time derivative of the true airspeed
        const auto TASdot_mps2 = (VTASxdot_mps2 * VTASx_mps + VTASydot_mps2 * VTASy_mps +
                                  VTASzdot_mps2 * VTASz_mps) / TAS_mps;


        // calculate the time derivatives of the
        // angles of attack and sideslip
        const auto [aoadot_radps,
                    aosdot_radps] = velocity_angles_time_derivatives(VTASxdot_mps2, VTASydot_mps2, VTASzdot_mps2,
                                                                     VTASx_mps, VTASy_mps, VTASz_mps);


        // calculate the angular velocity components of
        // the wind axes', with respect to Earth
        const auto [__ignore,
                    pqr_windREarth_radps] = velocity_axes_angular_velocities(p_radps, q_radps, r_radps,
                                                                             aoadot_radps, aosdot_radps,
                                                                             T_body2wind);
        (void)__ignore;


        // calculate the time derivatives of the air path angles
        const auto flight_path_angle_rad = tr7::deg2rad(flight_path_angle_deg);

        const auto [true_heading_angledot_radps,
                    flight_path_angledot_radps,
                    flight_roll_angledot_radps] = inverse_angular_kinematic_relationships(pqr_windREarth_radps[0], pqr_windREarth_radps[1], pqr_windREarth_radps[2],
                                                                                          tr7::deg2rad(true_heading_angle_deg),
                                                                                          flight_path_angle_rad,
                                                                                          tr7::deg2rad(flight_roll_angle_deg));


        // calculate the air turn radius
        const auto air_turn_radius_m = turn_radius(TAS_mps,
                                                   true_heading_angledot_radps,
                                                   flight_path_angledot_radps,
                                                   flight_path_angle_rad);


        // calculate the wind-axes load factors, i.e., the body-axes
        // ones but projected onto the wind-axes frame
        const auto [Nwx_g,
                    Nwy_g,
                    Nwz_g] = matrix3x3_times_vector3d(T_body2wind,
                                                      Nx_g, Ny_g, Nz_g);

        // return the model's outputs
        return outputs_type<T>{ dVTASdt_mps2[0], dVTASdt_mps2[1], dVTASdt_mps2[2],
                                VTASxdot_mps2, VTASydot_mps2, VTASzdot_mps2,
                                TASdot_mps2, tr7::rad2deg(aoadot_radps), tr7::rad2deg(aosdot_radps),
                                Nwx_g, Nwy_g, Nwz_g,
                                pqr_windREarth_radps[0], pqr_windREarth_radps[1], pqr_windREarth_radps[2],
                                tr7::rad2deg(true_heading_angledot_radps), tr7::rad2deg(flight_path_angledot_radps), tr7::rad2deg(flight_roll_angledot_radps),
                                air_turn_radius_m };

    }


}


#endif