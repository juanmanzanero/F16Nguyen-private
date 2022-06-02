#ifndef F16_NGUYEN_ALGEBRA3D
#define F16_NGUYEN_ALGEBRA3D
#pragma once


#include <array>

#include "tr7/tr7math.h"


//
// Defines handy functions that implement
// typical algebraic operations in 3d
// space, only employing scalar/std::array
// types.
//


namespace F16_Nguyen
{
    template<typename T>
    constexpr std::array<T, 3u> vector3d_normalize(const T &x, const T &y, const T &z)
    {
        //
        // Returns the normalized version of the input 3D vector.
        //

        const auto norm = tr7::AD_sqrt(x * x + y * y + z * z);

        const auto invnorm = norm > std::numeric_limits<T>::epsilon() ? T{ 1 } / norm : T{ 0 };

        return std::array<T, 3u>{ x * invnorm,
                                  y * invnorm,
                                  z * invnorm };

    }

    template<typename T>
    constexpr std::array<T, 3u> vector3d_normalize(const std::array<T, 3u> &v)
    {
        //
        // Returns the normalized version of the input 3D vector.
        //

        return vector3d_normalize(std::get<0>(v), std::get<1>(v), std::get<2>(v));

    }


    template<typename T>
    constexpr std::array<T, 3u> vector3d_cross_product(const T &lhsx, const T &lhsy, const T &lhsz,
                                                       const T &rhsx, const T &rhsy, const T &rhsz)
    {
        //
        // Returns the cross-product of the two input 3D vectors.
        //

        return std::array<T, 3u>{ lhsy * rhsz - lhsz * rhsy,
                                  lhsz * rhsx - lhsx * rhsz,
                                  lhsx * rhsy - lhsy * rhsx };

    }

    template<typename T>
    constexpr std::array<T, 3u> vector3d_cross_product(const T &lhsx, const T &lhsy, const T &lhsz,
                                                       const std::array<T, 3u> &rhs)
    {
        //
        // Returns the cross-product of the two input 3D vectors.
        //

        return vector3d_cross_product(lhsx, lhsy, lhsz,
                                      std::get<0>(rhs), std::get<1>(rhs), std::get<2>(rhs));

    }


    template<typename T>
    constexpr std::array<T, 9u> matrix3x3_transpose(const std::array<T, 9u> &M)
    {
        //
        // Returns the 3x3 matrix that results form transposing the
        // input 3x3 matrix. These matrices are stored as std::arrays
        // in column-major format.
        //

        return std::array<T, 9u>{ M[0],
                                  M[3],
                                  M[6],
                                  M[1],
                                  M[4],
                                  M[7],
                                  M[2],
                                  M[5],
                                  M[8] };

    }


    template<typename T>
    constexpr std::array<T, 9u> matrix3x3_times_matrix3x3(const std::array<T, 9u> &lhs,
                                                          const std::array<T, 9u> &rhs)
    {
        //
        // Returns the 3x3 matrix that results form multiplying
        // the two input 3x3 matrices. The 3 matrices are stored
        // as std::arrays in column-major format.
        //

        return std::array<T, 9u>{ lhs[0] * rhs[0] + lhs[3] * rhs[1] + lhs[6] * rhs[2],
                                  lhs[1] * rhs[0] + lhs[4] * rhs[1] + lhs[7] * rhs[2],
                                  lhs[2] * rhs[0] + lhs[5] * rhs[1] + lhs[8] * rhs[2],
                                  lhs[0] * rhs[3] + lhs[3] * rhs[4] + lhs[6] * rhs[5],
                                  lhs[1] * rhs[3] + lhs[4] * rhs[4] + lhs[7] * rhs[5],
                                  lhs[2] * rhs[3] + lhs[5] * rhs[4] + lhs[8] * rhs[5],
                                  lhs[0] * rhs[6] + lhs[3] * rhs[7] + lhs[6] * rhs[8],
                                  lhs[1] * rhs[6] + lhs[4] * rhs[7] + lhs[7] * rhs[8],
                                  lhs[2] * rhs[6] + lhs[5] * rhs[7] + lhs[8] * rhs[8] };

    }


    template<typename T>
    constexpr std::array<T, 3u> matrix3x3_times_vector3d(const std::array<T, 9u> &M,
                                                         const T &x, const T &y, const T &z)
    {
        //
        // Returns the 3d vector that results from multiplying
        // the input 3x3 matrix times the input 3-d vector. The
        // input matrix is an std::array given in column-major format.
        //

        return std::array<T, 3u>{ M[0] * x + M[3] * y + M[6] * z,
                                  M[1] * x + M[4] * y + M[7] * z,
                                  M[2] * x + M[5] * y + M[8] * z };

    }


    template<typename T>
    constexpr std::array<T, 3u> linsolve3x3(const std::array<T, 9u> &A,
                                            const T &bx, const T &by, const T &bz)
    {
        //
        // returns the solution "x" of the 3x3 linear
        // system "A * x = b". The input matrix should
        // come in col-major format.
        //

        const auto invdet = T{ 1 } / (A[0] * A[4] * A[8] +
                                      A[1] * A[5] * A[6] +
                                      A[3] * A[7] * A[2] -
                                      A[2] * A[6] * A[4] -
                                      A[5] * A[7] * A[0] -
                                      A[1] * A[3] * A[8]);

        return std::array<T, 3u>{ (bx * A[4] * A[8] +
                                   by * A[5] * A[6] +
                                   A[3] * A[7] * bz -
                                   bz * A[6] * A[4] -
                                   A[5] * A[7] * bx -
                                   by * A[3] * A[8]) * invdet,

                                  (A[0] * by * A[8] +
                                   A[1] * bz * A[6] +
                                   bx * A[7] * A[2] -
                                   A[2] * A[6] * by -
                                   bz * A[7] * A[0] -
                                   A[1] * bx * A[8]) * invdet,

                                  (A[0] * A[4] * bz +
                                   A[1] * A[5] * bx +
                                   A[3] * by * A[2] -
                                   A[2] * bx * A[4] -
                                   A[5] * by * A[0] -
                                   A[1] * A[3] * bz) * invdet };

    }


    template<typename T>
    constexpr std::array<T, 3u> rotmat2ea(const std::array<T, 9u> &M)
    {
        //
        // Returns an std::array holding the Euler angles in 'ZYX' sequence
        // "[yaw, pitch, roll]", in rad, from the input rotation matrix,
        // given as an std::array in column-major format.
        //

        // gymbal lock tolerance, calculated with
        // "gymbal_lock_pitch_tol_deg = 0.001;
        //  gymbal_lock_pitch_zone = cos(deg2rad(gymbal_lock_pitch_tol_deg));"
        constexpr_variable auto gymbal_lock_pitch_zone = T{ 0.999999999847691 };
        const auto sinp = -M[2];


        T pitch_rad;
        if (tr7::AD_abs(sinp) <= gymbal_lock_pitch_zone) {
            pitch_rad = tr7::AD_asin(sinp);
        }
        else if (sinp >= T{ 0 }) {
            pitch_rad = T{ 0.5 } * tr7::pi_T<T>;
        }
        else {
            pitch_rad = -T{ 0.5 } * tr7::pi_T<T>;
        }

        const auto yaw_rad = tr7::AD_atan2(M[1], M[0]);
        const auto roll_rad = tr7::AD_atan2(M[5], M[8]);

        return std::array<T, 3u>{ tr7::wrap_to_pi(yaw_rad),
                                  tr7::wrap_to_pi(pitch_rad),
                                  tr7::wrap_to_pi(roll_rad) };

    }


    template<typename T>
    constexpr std::array<T, 4u> quat_normalize(const T &qw, const T &qx, const T &qy, const T &qz)
    {
        //
        // Returns an std::array holding the normalized
        // version of input quaternion.
        //

        const auto norm = tr7::AD_sqrt(qw * qw + qx * qx + qy * qy + qz * qz);

        const auto invnorm = norm > std::numeric_limits<T>::epsilon() ? T{ 1 } / norm : T{ 0 };

        return std::array<T, 4u>{ qw * invnorm,
                                  qx * invnorm,
                                  qy * invnorm,
                                  qz * invnorm };

    }


    template<typename T>
    constexpr std::array<T, 9u> quat2rotmat(const T &qw, const T &qx, const T &qy, const T &qz)
    {
        //
        // Rotation matrix from quaternion. The rotation matrix (returned
        // as an std::array in column-major format) transforms from the
        // "rotated" axes to the "original" ones, like the input quaternion
        // does:
        //
        //    x_original      = rot_matrix * x_rotated
        //    [0; x_original] = quat_multiply(quat_multiply(q, [0; x_rotated]), quat_conjugate(q))
        //
        // NOTE: the input quaternion does not need to be normalized.
        //

        const auto xx = qx * qx;
        const auto yy = qy * qy;
        const auto zz = qz * qz;

        const auto two = T{ 2 } / (qw * qw + xx + yy + zz);

        const auto xx2 = two * xx;
        const auto yy2 = two * yy;
        const auto zz2 = two * zz;

        const auto wx2 = two * qw * qx;
        const auto wy2 = two * qw * qy;
        const auto wz2 = two * qw * qz;

        const auto xy2 = two * qx * qy;
        const auto xz2 = two * qx * qz;
        const auto yz2 = two * qy * qz;

        return std::array<T, 9u> { T{ 1 } - yy2 - zz2,
                                   xy2 + wz2,
                                   xz2 - wy2,
                                   xy2 - wz2,
                                   T{ 1 } - xx2 - zz2,
                                   yz2 + wx2,
                                   xz2 + wy2,
                                   yz2 - wx2,
                                   T{ 1 } - xx2 - yy2 };

    }


    template<typename T>
    constexpr std::array<T, 3u> quat2ea(const T &qw, const T &qx, const T &qy, const T &qz)
    {
        //
        // Returns an std::array holding the Euler angles in 'ZYX' sequence
        // "[yaw, pitch, roll]", in rad, from the input quaternion, which
        // doesn't need to be normalized.
        //

        // gymbal lock tolerance, calculated with
        // "gymbal_lock_pitch_tol_deg = 0.001;
        //  gymbal_lock_pitch_zone = cos(deg2rad(gymbal_lock_pitch_tol_deg));"
        constexpr_variable auto gymbal_lock_pitch_zone = T{ 0.999999999847691 };

        const auto ww = qw * qw;
        const auto xx = qx * qx;
        const auto yy = qy * qy;
        const auto zz = qz * qz;
        const auto two = T{ 2 } / (ww + xx + yy + zz);
        const auto sinp = two * (qw * qy - qz * qx);

        T pitch_rad;
        if (tr7::AD_abs(sinp) <= gymbal_lock_pitch_zone) {
            pitch_rad = tr7::AD_asin(sinp);
        }
        else if (sinp >= T{ 0 }) {
            pitch_rad = T{ 0.5 } * tr7::pi_T<T>;
        }
        else {
            pitch_rad = -T{ 0.5 } * tr7::pi_T<T>;
        }

        const auto yaw_rad = tr7::AD_atan2(T{ 2 } * (qw * qz + qx * qy), ww + xx - yy - zz);
        const auto roll_rad = tr7::AD_atan2(T{ 2 } * (qw * qx + qy * qz), ww - xx - yy + zz);

        return std::array<T, 3u>{ tr7::wrap_to_pi(yaw_rad),
                                  tr7::wrap_to_pi(pitch_rad),
                                  tr7::wrap_to_pi(roll_rad) };

    }


    template<typename T>
    constexpr std::array<T, 4u> quat_time_derivative(const T &qw, const T &qx, const T &qy, const T &qz,
                                                     const T &p_radpt, const T &q_radpt, const T &r_radpt)
    {
        //
        // Returns an std::array holding the time derivative of
        // the input quaterion ("qdot == dquat/dt") given the angular
        // velocity components "[p, q, r]" of a rotating frame (projected
        // in that same frame). The input quaternion MUST BE the
        // one that transforms FROM the rotating frame TO the
        // reference frame, and does not need to be normalized on
        // entry. The units of the result are [1 / time], with the
        // time being given by the units of the angular velocity inputs.
        //

        const auto [qnw, qnx, qny, qnz] = quat_normalize(qw, qx, qy, qz);

        return std::array<T, 4u>{ -T{ 0.5 } * (qnx * p_radpt + qny * q_radpt + qnz * r_radpt),
                                   T{ 0.5 } * (qnw * p_radpt + qny * r_radpt - qnz * q_radpt),
                                   T{ 0.5 } * (qnw * q_radpt - qnx * r_radpt + qnz * p_radpt),
                                   T{ 0.5 } * (qnw * r_radpt + qnx * q_radpt - qny * p_radpt) };

    }


    template<typename T>
    constexpr std::array<T, 3u> inverse_angular_kinematic_relationships(const T &p, const T &q, const T &r,
                                                                        const T &yaw_rad, const T &pitch_rad, const T &roll_rad)
    {
        //
        // Returns an std::array holding the time derivatives
        // of a body's Euler angles, i.e. "[yawdot, pitchdot, rolldot]",
        // from its body-axes angular velocity components "[p, q, r]"
        // and its Euler angles "[yaw_rad, pitch_rad, roll_rad]", in
        // the same units in which inputs "[p, q, r]" are given.
        //

        const auto cphi = tr7::AD_cos(roll_rad);
        const auto sphi = tr7::AD_sin(roll_rad);
        const auto rolldot = (q * sphi + r * cphi) / tr7::AD_cos(pitch_rad);

        return std::array<T, 3u>{ p + rolldot * tr7::AD_sin(pitch_rad),
                                  q * cphi - r * sphi,
                                  rolldot };

    }


    template<typename T>
    constexpr std::array<T, 2> velocity_angles(const T &u, const T &v, const T &w)
    {
        //
        // Returns an std::array holding a body's "[alpha, beta]" [rad]
        // angles given its body-axes velocity components "[u, v, w]".
        // For example, if this velocity vector is the body's VTAS (true
        // airspeed), "alpha" and "beta" will be the body's aerodynamic
        // angles of attack and sideslip.
        //

        const auto lon2 = u * u + w * w;
        //if (lon2 + v * v > std::numeric_limits<T>::epsilon()) {
            const auto sgn_u = u >= T{ 0 } ? T{ 1 } : T{ -1 };

            const auto alpha_rad = tr7::AD_atan2(sgn_u * w, tr7::AD_abs(u));
            const auto beta_rad = tr7::wrap_to_pi(tr7::AD_atan2(v, sgn_u * tr7::AD_sqrt(lon2)));

            return std::array<T, 2>{ alpha_rad, beta_rad };

        //}
        //else {
        //    // zero velocity convention
        //    return std::array<T, 2>{ 0, 0 };
        //
        //}

    }


    template<typename T>
    constexpr std::array<T, 9u> velocity_angles2rotmat(const T &alpha_rad, const T &beta_rad)
    {
        //
        // Returns the "T_body2velaxes" rotation matrix, i.e.,
        // "x_velocityaxes = T_body2velaxes * x_bodyaxes", from
        // the input velocity angles (see function "velocity_angles").
        // For example, if the input velocity angles are the body's
        // aerodynamic angles of attack and sideslip, "T_body2velocity"
        // will be the famous "T_body2wind" rotation matrix that
        // transforms from the body-axes frame to the wind-axes one.
        //
        // NOTE: the output matrix is returned as an std::array
        // in column-major format.
        //

        const auto calpha = tr7::AD_cos(alpha_rad);
        const auto salpha = tr7::AD_sin(alpha_rad);

        const auto cbeta = tr7::AD_cos(beta_rad);
        const auto sbeta = tr7::AD_sin(beta_rad);

        return std::array<T, 9u>{ calpha * cbeta,
                                  -calpha * sbeta,
                                  -salpha,
                                  sbeta,
                                  cbeta,
                                  T{ 0 },
                                  salpha * cbeta,
                                  -salpha * sbeta,
                                  calpha };

    }


    template<typename T>
    constexpr std::array<T, 2u> velocity_angles_time_derivatives(const T &udot, const T &vdot, const T &wdot,
                                                                 const T &u, const T &v, const T &w)
    {
        //
        // Returns an std::array holding "[alphadot, betadot]", i.e., the
        // time derivatives of a body's velocity angles, in [rad / time], from
        // the input pseudo-acceleration components "[udot, vdot, wdot]"
        // and velocity components "[u, v, w]", both projected in body axes.
        //
        // For example, if the input velocity and pseudo-acceleration components
        // come from the body's VTAS (true airspeed vector), "alphadot" and "betadot"
        // will be the time derivatives of the body's aerodynamic angles of attack
        // and sideslip.
        //
        // NOTE: what we call "pseudo-acceleration" is the explicit
        // time-derivative in the body-axes frame of a body's velocity, i.e.
        //
        // "acceleration = pseudo_acceleration + cross(pqr, velocity),"
        //
        // all in the body-axes frame.
        //

        const auto sgn_u = u >= T{ 0 } ? T{ 1 } : T{ -1 };
        const auto lon2 = u * u + w * w;
        const auto den = lon2 > std::numeric_limits<T>::epsilon() ? lon2 : std::numeric_limits<T>::epsilon();

        const auto alphadot_radpt = (wdot * u - udot * w) / den;
        const auto betadot_radpt = sgn_u * (vdot * lon2 - v * (udot * u + wdot * w)) /
                                   ((den + v * v) * tr7::AD_sqrt(den));

        return std::array<T, 2u>{ alphadot_radpt, betadot_radpt };

    }


    template<typename T>
    constexpr std::pair<std::array<T, 3u>,
                        std::array<T, 3u> > velocity_axes_angular_velocities(const T &p, const T &q, const T &r,
                                                                             const T &alphadot, const T &betadot,
                                                                             const std::array<T, 9u> &T_body2velaxes)
    {
        //
        // Returns an std::pair holding two angular velocity
        // vectors, "[pqr_bodyRvelaxes, pqr_velaxesRref]",
        // in the same units as inputs "p, q, r, alphadot,
        // betadot".
        //
        // "pqr_bodyRvelaxes" is the angular velocity of the
        // body-axes frame relative to the velocity-axes
        // one, projected on body-axes. "pqr_velaxesRREF"
        // is the angular velocity of the velocity-axes relative
        // to the reference frame (e.g. the Earth-axes), expressed
        // with respect to the velocity-axes frame.
        //
        // Inputs "[p, q, r]" are the body-axes components of the
        // body's angular velocity with respect to the reference
        // frame, "[alphadot, betadot]" are the time derivatives
        // of the body's velocity angles (see function
        // "velocity_angles_time_derivatives"), and "T_body2velaxes"
        // is the col-major rotation matrix formed with these
        // velocity angles (see function "velocity_angles2rotmat").
        //

        const auto salpha = -T_body2velaxes[2];
        const auto &calpha = T_body2velaxes[8];

        // this one is just "-betadot * z_wind + alphadot * y_body"
        const auto pqr_bodyRvelaxes = std::array<T, 3u> { betadot * salpha,
                                                          alphadot,
                                                          -betadot * calpha };

        // this one has the [pw, qw, rw] components
        const auto pqr_velaxesREarth = matrix3x3_times_vector3d(T_body2velaxes,
                                                                p - pqr_bodyRvelaxes[0],
                                                                q - pqr_bodyRvelaxes[1],
                                                                r - pqr_bodyRvelaxes[2]);

        return std::make_pair(pqr_bodyRvelaxes, pqr_velaxesREarth);

    }


    template<typename T>
    constexpr std::array<T, 3u> pseudoacceleration2acceleration(const T &udot, const T &vdot, const T &wdot,
                                                                const T &u, const T &v, const T &w,
                                                                const T &p_radpt, const T &q_radpt, const T &r_radpt)
    {
        //
        // Returns an std::array holding the body-axes
        // acceleration components "[ax, ay, az]"
        // from the input body-axes pseudo-acceleration,
        // velocity and angular speed components. The units
        // of the result are those of the inputs
        // "[udot, vdot, wdot]", and velocities
        // "[u, v, w]", "[p, q, r]" should be compatible
        // in time and length with the former.
        //
        // NOTE: what we call "pseudo-acceleration" is
        // the explicit time-derivative in the body-axes
        // frame of a body's velocity, i.e.
        //
        // "acceleration = pseudo_acceleration + cross(pqr, velocity),"
        //
        // all in the body-axes frame.
        //

        return std::array<T, 3u>{ udot + q_radpt * w - r_radpt * v,
                                  vdot + r_radpt * u - p_radpt * w,
                                  wdot + p_radpt * v - q_radpt * u };

    }


    template<typename T>
    constexpr std::array<T, 3u> acceleration2pseudoacceleration(const T &ax, const T &ay, const T &az,
                                                                const T &u, const T &v, const T &w,
                                                                const T &p_radpt, const T &q_radpt, const T &r_radpt)
    {
        //
        // Returns an std::array holding the body-axes
        // pseudo-acceleration components "[udot, vdot, wdot]"
        // from the input body-axes acceleration, velocity
        // and angular speed components. The units of the
        // results are those of the inputs "[ax, ay, az]",
        // and velocities "[u, v, w]", "[p, q, r]" should
        // be compatible in time and length with the former.
        //
        // NOTE: what we call "pseudo-acceleration" is
        // the explicit time-derivative in the body-axes
        // frame of a body's velocity, i.e.
        //
        // "acceleration = pseudo_acceleration + cross(pqr, velocity),"
        //
        // all in the body-axes frame.
        //

        return pseudoacceleration2acceleration(ax, ay, az,
                                               u, v, w,
                                               -p_radpt, -q_radpt, -r_radpt);

    }


    template<typename T>
    constexpr T turn_radius(const T &speed, const T &horizontal_heading_angledot_radps,
                            const T &vertical_path_angledot_radps, const T &vertical_path_angle_rad)
    {
        //
        // Returns the turn radius of a trajectory, with respect
        // to a reference frame (e.g., the Earth axes one) given
        // the speed and the (time derivatives of the) angles that
        // determine the velocity vector in the reference frame.
        // The result is given in the same length units that the
        // input "speed" has. For example, if the input speed is
        // the true airspeed, the "horizontal_headingdot_radps" has
        // to be "chidot" (time derivative of the "true heading"), 
        // the "vertical_path_angle" is "gamma" (a.k.a "flight path
        // angle"), and the "vertical_path_angledot_radps" its time
        // derivative "gammadot", for the final result to be the
        // "air turn radius".
        //

        const auto chidot_times_cgamma = horizontal_heading_angledot_radps *
                                         tr7::AD_cos(vertical_path_angle_rad);

        return speed / tr7::AD_sqrt(chidot_times_cgamma * chidot_times_cgamma +
                                    vertical_path_angledot_radps * vertical_path_angledot_radps);

    }


}


#endif
