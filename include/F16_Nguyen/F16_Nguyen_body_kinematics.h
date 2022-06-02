#ifndef F16_NGUYEN_F16_NGUYEN_BODY_KINEMATICS_H
#define F16_NGUYEN_F16_NGUYEN_BODY_KINEMATICS_H
#pragma once


#include "F16_Nguyen/modelglue.h"
#include "F16_Nguyen/algebra3d.h"


//
// Defines the "F16_Nguyen_body_kinematics" class,
// that calculates kinematic variables of the
// body from its state.
//


namespace F16_Nguyen
{
    class F16_Nguyen_body_kinematics
    {
        //
        // Calculates several body-frame-related
        // kinematic variables.
        //

    public:

        DECLARE_MODEL_OUTPUTS(VGSx_Earthaxes_mps, VGSy_Earthaxes_mps, VGSz_Earthaxes_mps,
                              DOWNx_bodyaxes, DOWNy_bodyaxes, DOWNz_bodyaxes,
                              yawdot_degps, pitchdot_degps, rolldot_degps,
                              yaw_deg, pitch_deg, roll_deg,
                              GS_mps, ground_aoa_deg, ground_aos_deg,
                              ground_true_heading_angle_deg, ground_flight_path_angle_deg, ground_flight_roll_angle_deg)


        template<typename T>
        static constexpr outputs_type<T> outputs(const T &VGSx_mps, const T &VGSy_mps, const T &VGSz_mps,
                                                 const T &p_radps, const T &q_radps, const T &r_radps,
                                                 const T &qw_body2Earth, const T &qx_body2Earth, const T &qy_body2Earth, const T &qz_body2Earth);

    };


    //
    // "F16_Nguyen_body_kinematics" impl.
    //

    template<typename T>
    constexpr auto F16_Nguyen_body_kinematics::outputs(const T &VGSx_mps, const T &VGSy_mps, const T &VGSz_mps,
                                                       const T &p_radps, const T &q_radps, const T &r_radps,
                                                       const T &qw_body2Earth, const T &qx_body2Earth, const T &qy_body2Earth, const T &qz_body2Earth) -> outputs_type<T>
    {
        //
        // Evaluates the model's output equation, i.e.,
        // "outputs = G(state, inputs)", and returns the
        // model's "y".
        //

        // calculate the body-to-Earth rotation matrix
        const auto T_body2Earth = quat2rotmat(qw_body2Earth, qx_body2Earth, qy_body2Earth, qz_body2Earth);


        // calculate the groundspeed vector and norm, in Earthaxes
        const auto [VGSx_Earthaxes_mps,
                    VGSy_Earthaxes_mps,
                    VGSz_Earthaxes_mps] = matrix3x3_times_vector3d(T_body2Earth,
                                                                   VGSx_mps, VGSy_mps, VGSz_mps);


        // calculate the "DOWN" vector in bodyaxes
        // ("DOWN" == "z_Earth"), it determines the
        // direction of the gravity
        const auto &DOWNx_bodyaxes = T_body2Earth[2];
        const auto &DOWNy_bodyaxes = T_body2Earth[5];
        const auto &DOWNz_bodyaxes = T_body2Earth[8];


        // calculate the Euler angles
        const auto [yaw_rad,
                    pitch_rad,
                    roll_rad] = quat2ea(qw_body2Earth, qx_body2Earth, qy_body2Earth, qz_body2Earth);


        // calculate the Euler angle rates
        const auto [yawdot_radps,
                    pitchdot_radps,
                    rolldot_radps] = inverse_angular_kinematic_relationships(p_radps, q_radps, r_radps,
                                                                             yaw_rad, pitch_rad, roll_rad);


        // calculate the groundspeed
        const auto GS_mps = tr7::AD_sqrt(VGSx_mps * VGSx_mps +
                                         VGSy_mps * VGSy_mps +
                                         VGSz_mps * VGSz_mps);


        // calculate the ground angles of attack, sideslip, and path angles
        const auto [ground_aoa_rad,
                    ground_aos_rad] = velocity_angles(VGSx_mps, VGSy_mps, VGSz_mps);

        const auto T_body2groundvel = velocity_angles2rotmat(ground_aoa_rad, ground_aos_rad);
        const auto T_groundvel2Earth = matrix3x3_times_matrix3x3(T_body2Earth, matrix3x3_transpose(T_body2groundvel));

        const auto [ground_true_heading_angle_rad,
                    ground_flight_path_angle_rad,
                    ground_flight_roll_angle_rad] = rotmat2ea(T_groundvel2Earth);


        // return the model's outputs
        return outputs_type<T>{ VGSx_Earthaxes_mps, VGSy_Earthaxes_mps, VGSz_Earthaxes_mps,
                                DOWNx_bodyaxes, DOWNy_bodyaxes, DOWNz_bodyaxes,
                                tr7::rad2deg(yawdot_radps), tr7::rad2deg(pitchdot_radps), tr7::rad2deg(rolldot_radps),
                                tr7::rad2deg(yaw_rad), tr7::rad2deg(pitch_rad), tr7::rad2deg(roll_rad),
                                GS_mps, tr7::rad2deg(ground_aoa_rad), tr7::rad2deg(ground_aos_rad),
                                tr7::rad2deg(ground_true_heading_angle_rad), tr7::rad2deg(ground_flight_path_angle_rad), tr7::rad2deg(ground_flight_roll_angle_rad) };

    }


}


#endif