#ifndef F16_NGUYEN_F16_NGUYEN_FLAT_EARTH_BODY_DYNAMICS_H
#define F16_NGUYEN_F16_NGUYEN_FLAT_EARTH_BODY_DYNAMICS_H
#pragma once


#include "F16_Nguyen/modelglue.h"
#include "F16_Nguyen/contiguous_ranges.h"
#include "F16_Nguyen/F16_Nguyen_body_kinematics.h"


//
// Defines the "F16_Nguyen_flat_Earth_body_dynamics" class,
// that implements the movement equations for the simulation
// model explained in reference "Simulator study of stall/post-
// stall characteristics of a fighter airplane with relaxed
// longitudinal static stability", L.T.Nguyen et al.,
// NASA-TP-1538, 1979; employing a flat-Earth approximation.
//


namespace F16_Nguyen
{
    class F16_Nguyen_flat_Earth_body_dynamics
    {
        //
        // A class that implements the movement
        // equations for Nguyen's F16 model, using
        // the flat-Earth approximation.
        //

    public:

        using data_type = double;


        struct default_parameters
        {
            static constexpr data_type mass_kg = 91188 / 9.81;
            static constexpr data_type Ixx_kgm2 = 12875.;
            static constexpr data_type Ixy_kgm2 = 0.;
            static constexpr data_type Ixz_kgm2 = 1331.;
            static constexpr data_type Iyy_kgm2 = 75674.;
            static constexpr data_type Iyz_kgm2 = 0.;
            static constexpr data_type Izz_kgm2 = 85552.;

            static constexpr data_type thrust_point_of_applicationx_bodyaxes_m = -6.; // guesswork...
            static constexpr data_type thrust_point_of_applicationy_bodyaxes_m = 0.;
            static constexpr data_type thrust_point_of_applicationz_bodyaxes_m = 0.;

        };


        DECLARE_MODEL_STATES(u_mps, v_mps, w_mps,
                             p_radps, q_radps, r_radps,
                             qw_body2Earth, qx_body2Earth, qy_body2Earth, qz_body2Earth,
                             xEarth_m, yEarth_m, zEarth_m)

        DECLARE_MODEL_DERIVED_OUTPUTS(F16_Nguyen_body_kinematics, Ftotx_N, Ftoty_N, Ftotz_N,
                                                                  Mtotx_N, Mtoty_N, Mtotz_N,
                                                                  Nx_g, Ny_g, Nz_g,
                                                                  ax_mps2, ay_mps2, az_mps2,
                                                                  udot_mps2, vdot_mps2, wdot_mps2,
                                                                  pdot_radps2, qdot_radps2, rdot_radps2,
                                                                  GSdot_mps, ground_aoadot_degps, ground_aosdot_degps,
                                                                  ground_true_heading_angledot_degps, ground_flight_path_angledot_degps, ground_flight_roll_angledot_degps,
                                                                  ground_turn_radius_m)


        template<typename StatesType, typename T>
        static constexpr outputs_type<T> outputs(const StatesType &x,
                                                 const T &Faerox_N, const T &Faeroy_N, const T &Faeroz_N,
                                                 const T &Maerox_Nm, const T &Maeroy_Nm, const T &Maeroz_Nm,
                                                 const T &engine_thrust_N, const T &engine_angular_momentum_kgm2ps,
                                                 const T &thrustvectoring_longitude_angle_deg, const T &thrustvectoring_latitude_angle_deg,
                                                 const T &mass_kg = T{ default_parameters::mass_kg },
                                                 const T &Ixx_kgm2 = T{ default_parameters::Ixx_kgm2 },
                                                 const T &Ixy_kgm2 = T{ default_parameters::Ixy_kgm2 },
                                                 const T &Ixz_kgm2 = T{ default_parameters::Ixz_kgm2 },
                                                 const T &Iyy_kgm2 = T{ default_parameters::Iyy_kgm2 },
                                                 const T &Iyz_kgm2 = T{ default_parameters::Iyz_kgm2 },
                                                 const T &Izz_kgm2 = T{ default_parameters::Izz_kgm2 },
                                                 const T &thrust_point_of_applicationx_bodyaxes_m = T{ default_parameters::thrust_point_of_applicationx_bodyaxes_m },
                                                 const T &thrust_point_of_applicationy_bodyaxes_m = T{ default_parameters::thrust_point_of_applicationy_bodyaxes_m },
                                                 const T &thrust_point_of_applicationz_bodyaxes_m = T{ default_parameters::thrust_point_of_applicationz_bodyaxes_m });

        template<typename StatesType, typename OutputsType,
                 typename T = tr7::typeof_bracket_operator<StatesType> >
        static constexpr states_type<T> derivatives(const StatesType &x,
                                                    const OutputsType &y);

        template<typename StatesType,
                 typename T = tr7::typeof_bracket_operator<StatesType> >
        static constexpr states_type<T> state_limiters(const StatesType &x);

    };


    //
    // "F16_Nguyen_flat_Earth_body_dynamics" impl.
    //

    template<typename StatesType, typename T>
    constexpr auto F16_Nguyen_flat_Earth_body_dynamics::outputs(const StatesType &x,
                                                                const T &Faerox_N, const T &Faeroy_N, const T &Faeroz_N,
                                                                const T &Maerox_Nm, const T &Maeroy_Nm, const T &Maeroz_Nm,
                                                                const T &engine_thrust_N, const T &engine_angular_momentum_kgm2ps,
                                                                const T &thrustvectoring_longitude_angle_deg, const T &thrustvectoring_latitude_angle_deg,
                                                                const T &mass_kg,
                                                                const T &Ixx_kgm2,
                                                                const T &Ixy_kgm2,
                                                                const T &Ixz_kgm2,
                                                                const T &Iyy_kgm2,
                                                                const T &Iyz_kgm2,
                                                                const T &Izz_kgm2,
                                                                const T &thrust_point_of_applicationx_bodyaxes_m,
                                                                const T &thrust_point_of_applicationy_bodyaxes_m,
                                                                const T &thrust_point_of_applicationz_bodyaxes_m) -> outputs_type<T>
    {
        //
        // Evaluates the model's output equation, i.e.,
        // "outputs = G(state, inputs)", and returns the
        // model's "y".
        //

        ensure_minimum_size<_num_states>(x);
        outputs_type<T> y{};


        // start by calculating the base outputs
        std::copy_n(F16_Nguyen_body_kinematics::outputs(x[state_names::u_mps], x[state_names::v_mps], x[state_names::w_mps],
                                                        x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps],
                                                        x[state_names::qw_body2Earth], x[state_names::qx_body2Earth], x[state_names::qy_body2Earth], x[state_names::qz_body2Earth]).cbegin(),
                    F16_Nguyen_body_kinematics::_num_outputs, std::begin(y));


        // calculate the "engine line", along which we'll apply its thrust
        // and angular momentum. It is a unit vector defined by its
        // [longitude, latitude] angles with respect to the body axes frame
        const auto thrustvectoring_longitude_angle_rad = tr7::deg2rad(thrustvectoring_longitude_angle_deg);
        const auto thrustvectoring_latitude_angle_rad = tr7::deg2rad(thrustvectoring_latitude_angle_deg);
        const auto cenginelat = tr7::AD_cos(thrustvectoring_latitude_angle_rad);

        const auto engine_line_x = cenginelat * tr7::AD_cos(thrustvectoring_longitude_angle_rad);
        const auto engine_line_y = cenginelat * tr7::AD_sin(thrustvectoring_longitude_angle_rad);
        const auto engine_line_z = tr7::AD_sin(thrustvectoring_latitude_angle_rad);


        // calculate the engine forces
        const auto Fenginex_N = engine_thrust_N * engine_line_x;
        const auto Fenginey_N = engine_thrust_N * engine_line_y;
        const auto Fenginez_N = engine_thrust_N * engine_line_z;


        // calculate the engine moments: cross(angular_momentum * engine_line, pqr) +
        // cross(thrust_point_of_application, thrust * engine_line)
        const auto [Mengine_angular_momentumx_Nm,
                    Mengine_angular_momentumy_Nm,
                    Mengine_angular_momentumz_Nm] = vector3d_cross_product(engine_angular_momentum_kgm2ps * engine_line_x,
                                                                           engine_angular_momentum_kgm2ps * engine_line_y,
                                                                           engine_angular_momentum_kgm2ps * engine_line_z,
                                                                           x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps]);

        const auto [Mengine_thrustx_Nm,
                    Mengine_thrusty_Nm,
                    Mengine_thrustz_Nm] = vector3d_cross_product(thrust_point_of_applicationx_bodyaxes_m,
                                                                 thrust_point_of_applicationy_bodyaxes_m,
                                                                 thrust_point_of_applicationz_bodyaxes_m,
                                                                 Fenginex_N, Fenginey_N, Fenginez_N);


        // calculate the total forces & moments (note: we could add the "engine_angular_momentumdot" term...)
        y[output_names::Ftotx_N] = Faerox_N + Fenginex_N;
        y[output_names::Ftoty_N] = Faeroy_N + Fenginey_N;
        y[output_names::Ftotz_N] = Faeroz_N + Fenginez_N;

        y[output_names::Mtotx_N] = Maerox_Nm + Mengine_angular_momentumx_Nm + Mengine_thrustx_Nm;
        y[output_names::Mtoty_N] = Maeroy_Nm + Mengine_angular_momentumy_Nm + Mengine_thrusty_Nm;
        y[output_names::Mtotz_N] = Maeroz_Nm + Mengine_angular_momentumz_Nm + Mengine_thrustz_Nm;


        // calculate the load factors & acceleration
        constexpr_variable auto g0_mps2 = T{ 9.80665 };

        const auto weight_N = mass_kg * g0_mps2;
        y[output_names::Nx_g] = -y[output_names::Ftotx_N] / weight_N;
        y[output_names::Ny_g] = -y[output_names::Ftoty_N] / weight_N;
        y[output_names::Nz_g] = -y[output_names::Ftotz_N] / weight_N;

        y[output_names::ax_mps2] = g0_mps2 * (y[output_names::DOWNx_bodyaxes] - y[output_names::Nx_g]);
        y[output_names::ay_mps2] = g0_mps2 * (y[output_names::DOWNy_bodyaxes] - y[output_names::Ny_g]);
        y[output_names::az_mps2] = g0_mps2 * (y[output_names::DOWNz_bodyaxes] - y[output_names::Nz_g]);


        // calculate the linear pseudo-acceleration
        const auto [udot_mps2,
                    vdot_mps2,
                    wdot_mps2] = acceleration2pseudoacceleration(y[output_names::ax_mps2], y[output_names::ay_mps2], y[output_names::az_mps2],
                                                                 x[state_names::u_mps], x[state_names::v_mps], x[state_names::w_mps],
                                                                 x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps]);

        y[output_names::udot_mps2] = udot_mps2;
        y[output_names::vdot_mps2] = vdot_mps2;
        y[output_names::wdot_mps2] = wdot_mps2;


        // calculate the angular acceleration
        const auto inertia_matrix_kgm2 = std::array<T, 9u>{ Ixx_kgm2, Ixy_kgm2, Ixz_kgm2,
                                                            Ixy_kgm2, Iyy_kgm2, Iyz_kgm2,
                                                            Ixz_kgm2, Iyz_kgm2, Izz_kgm2 };

        const auto [wIw_x,
                    wIw_y,
                    wIw_z] = vector3d_cross_product(x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps],
                                                   matrix3x3_times_vector3d(inertia_matrix_kgm2,
                                                                            x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps]));

        const auto [pdot_radps2,
                    qdot_radps2,
                    rdot_radps2] = linsolve3x3(inertia_matrix_kgm2,
                                               y[output_names::Mtotx_N] - wIw_x,
                                               y[output_names::Mtoty_N] - wIw_y,
                                               y[output_names::Mtotz_N] - wIw_z);

        y[output_names::pdot_radps2] = pdot_radps2;
        y[output_names::qdot_radps2] = qdot_radps2;
        y[output_names::rdot_radps2] = rdot_radps2;


        // calculate the time derivative of the groundspeed
        y[output_names::GSdot_mps] = (y[output_names::udot_mps2] * x[state_names::u_mps] +
                                      y[output_names::vdot_mps2] * x[state_names::v_mps] +
                                      y[output_names::wdot_mps2] * x[state_names::w_mps]) / y[output_names::GS_mps];


        // calculate the derivatives of the ground angles
        // of attack, sideslip, and path angles
        const auto [ground_aoadot_radps,
                    ground_aosdot_radps] = velocity_angles_time_derivatives(y[output_names::udot_mps2], y[output_names::vdot_mps2], y[output_names::wdot_mps2],
                                                                            x[state_names::u_mps], x[state_names::v_mps], x[state_names::w_mps]);

        y[output_names::ground_aoadot_degps] = tr7::rad2deg(ground_aoadot_radps);
        y[output_names::ground_aosdot_degps] = tr7::rad2deg(ground_aosdot_radps);

        const auto [__ignore,
                    pqr_GSREarth_radps] = velocity_axes_angular_velocities(x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps],
                                                                           ground_aoadot_radps, ground_aosdot_radps,
                                                                           velocity_angles2rotmat(tr7::deg2rad(y[output_names::ground_aoa_deg]),
                                                                                                  tr7::deg2rad(y[output_names::ground_aos_deg])));
        (void)__ignore;

        const auto [ground_true_heading_angledot_radps,
                    ground_flight_path_angledot_radps,
                    ground_flight_roll_angledot_radps] = inverse_angular_kinematic_relationships(pqr_GSREarth_radps[0], pqr_GSREarth_radps[1], pqr_GSREarth_radps[2],
                                                                                                 tr7::deg2rad(y[output_names::ground_true_heading_angle_deg]),
                                                                                                 tr7::deg2rad(y[output_names::ground_flight_path_angle_deg]),
                                                                                                 tr7::deg2rad(y[output_names::ground_flight_roll_angle_deg]));

        y[output_names::ground_true_heading_angledot_degps] = tr7::rad2deg(ground_true_heading_angledot_radps);
        y[output_names::ground_flight_path_angledot_degps] = tr7::rad2deg(ground_flight_path_angledot_radps);
        y[output_names::ground_flight_roll_angledot_degps] = tr7::rad2deg(ground_flight_roll_angledot_radps);


        // calculate the ground turn radius
        y[output_names::ground_turn_radius_m] = turn_radius(y[output_names::GS_mps],
                                                            ground_true_heading_angledot_radps,
                                                            ground_flight_path_angledot_radps,
                                                            tr7::deg2rad(y[output_names::ground_flight_path_angle_deg]));


        // return the model's outputs
        return y;

    }


    template<typename StatesType, typename OutputsType,
             typename T>
    constexpr auto F16_Nguyen_flat_Earth_body_dynamics::derivatives(const StatesType &x,
                                                                    const OutputsType &y) -> states_type<T>
    {
        //
        // Evaluates the model's state equation, i.e.,
        // "statedot = F(state, outputs, ...)", and returns
        // the model's "xdot".
        //

        ensure_minimum_size<_num_states>(x);
        ensure_minimum_size<_num_outputs>(y);
        states_type<T> xdot{};


        // linear pseudo-accelerations
        xdot[state_names::u_mps] = y[output_names::udot_mps2];
        xdot[state_names::v_mps] = y[output_names::vdot_mps2];
        xdot[state_names::w_mps] = y[output_names::wdot_mps2];


        // angular accelerations
        xdot[state_names::p_radps] = y[output_names::pdot_radps2];
        xdot[state_names::q_radps] = y[output_names::qdot_radps2];
        xdot[state_names::r_radps] = y[output_names::rdot_radps2];


        // quaterniondot
        const auto [qdotw_body2Earth,
                    qdotx_body2Earth,
                    qdoty_body2Earth,
                    qdotz_body2Earth] = quat_time_derivative(x[state_names::qw_body2Earth], x[state_names::qx_body2Earth], x[state_names::qy_body2Earth], x[state_names::qz_body2Earth],
                                                             x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps]);

        xdot[state_names::qw_body2Earth] = qdotw_body2Earth;
        xdot[state_names::qx_body2Earth] = qdotx_body2Earth;
        xdot[state_names::qy_body2Earth] = qdoty_body2Earth;
        xdot[state_names::qz_body2Earth] = qdotz_body2Earth;


        // positiondot
        xdot[state_names::xEarth_m] = y[output_names::VGSx_Earthaxes_mps];
        xdot[state_names::yEarth_m] = y[output_names::VGSy_Earthaxes_mps];
        xdot[state_names::zEarth_m] = y[output_names::VGSz_Earthaxes_mps];


        // return the model's statedots
        return xdot;

    }


    template<typename StatesType,
             typename T>
    constexpr auto F16_Nguyen_flat_Earth_body_dynamics::state_limiters(const StatesType &x) -> states_type<T>
    {
        //
        // Evaluates the model's state limiter equation, i.e.,
        // "x = SL(x, ...)", and returns the model's limited "x".
        //

        ensure_minimum_size<_num_states>(x);
        states_type<T> xlim{};


        // copy the input state
        tr7::tuple_index_for(xlim, [&x](auto &xlim_i, auto i) { xlim_i = x[i]; });


        // normalize the quaternion
        const auto[qnw_body2Earth,
                   qnx_body2Earth,
                   qny_body2Earth,
                   qnz_body2Earth] = quat_normalize(x[state_names::qw_body2Earth],
                                                    x[state_names::qx_body2Earth],
                                                    x[state_names::qy_body2Earth],
                                                    x[state_names::qz_body2Earth]);

        xlim[state_names::qw_body2Earth] = qnw_body2Earth;
        xlim[state_names::qx_body2Earth] = qnx_body2Earth;
        xlim[state_names::qy_body2Earth] = qny_body2Earth;
        xlim[state_names::qz_body2Earth] = qnz_body2Earth; 


        // return the limited state
        return xlim;

    }


}


#endif