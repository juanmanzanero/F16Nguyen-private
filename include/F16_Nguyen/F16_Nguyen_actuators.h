#ifndef F16_NGUYEN_F16_NGUYEN_ACTUATORS_H
#define F16_NGUYEN_F16_NGUYEN_ACTUATORS_H
#pragma once


#include "tr7/tr7tuple.h"
#include "tr7/tr7math.h"

#include "F16_Nguyen/modelglue.h"
#include "F16_Nguyen/contiguous_ranges.h"


//
// Defines the "F16_Nguyen_actuators" class, implementing
// the actuator model explained in reference "Simulator study of
// stall/post-stall characteristics of a fighter airplane with
// relaxed longitudinal static stability", L.T.Nguyen et al.,
// NASA-TP-1538, 1979.
//


namespace F16_Nguyen
{
    class F16_Nguyen_actuators
    {
        //
        // A class to evaluate Nguyen's F16 actuator model.
        // Its inputs are the actuator position demands, and its
        // states represent the actual surface deflections. All
        // actuators are modelled via first-order low pass
        // filters, with different rate and position limits.
        //

    public:

        using data_type = double;


        struct default_parameters
        {
            static constexpr data_type tau_stabilator_s = 0.0495;
            static constexpr data_type tau_leading_edge_flap_s = 0.136;
            static constexpr data_type tau_speedbreak_s = tau_stabilator_s;
            static constexpr data_type tau_aileron_s = tau_stabilator_s;
            static constexpr data_type tau_rudder_s = tau_stabilator_s;
            static constexpr data_type tau_thrustvectoring_longitude_angle_s = tau_leading_edge_flap_s;
            static constexpr data_type tau_thrustvectoring_latitude_angle_s = tau_thrustvectoring_longitude_angle_s;

            static constexpr data_type maxabs_stabilator_rate_degps = 60.;
            static constexpr data_type maxabs_leading_edge_flap_rate_degps = 25.;
            static constexpr data_type maxabs_speedbreak_rate_degps = 120.;
            static constexpr data_type maxabs_aileron_rate_degps = 80.;
            static constexpr data_type maxabs_rudder_rate_degps = maxabs_speedbreak_rate_degps;
            static constexpr data_type maxabs_thrustvectoring_longitude_angle_rate_degps = maxabs_leading_edge_flap_rate_degps;
            static constexpr data_type maxabs_thrustvectoring_latitude_angle_rate_degps = maxabs_thrustvectoring_longitude_angle_rate_degps;

            static constexpr data_type max_stabilator_deg = 25.;
            static constexpr data_type max_leading_edge_flap_deg = 25.;
            static constexpr data_type max_speedbreak_deg = 60.;
            static constexpr data_type max_aileron_deg = 21.5;
            static constexpr data_type max_rudder_deg = 30.;
            static constexpr data_type max_thrustvectoring_longitude_angle_deg = max_stabilator_deg;
            static constexpr data_type max_thrustvectoring_latitude_angle_deg = max_thrustvectoring_longitude_angle_deg;

            static constexpr data_type min_stabilator_deg = -max_stabilator_deg;
            static constexpr data_type min_leading_edge_flap_deg = 0.;
            static constexpr data_type min_speedbreak_deg = 0.;
            static constexpr data_type min_aileron_deg = -max_aileron_deg;
            static constexpr data_type min_rudder_deg = -max_rudder_deg;
            static constexpr data_type min_thrustvectoring_longitude_angle_deg = -max_thrustvectoring_longitude_angle_deg;
            static constexpr data_type min_thrustvectoring_latitude_angle_deg = -max_thrustvectoring_latitude_angle_deg;

        };


        DECLARE_MODEL_STATES(dh_deg, dlef_deg, dsb_deg,
                             da_deg, dr_deg,
                             thrustvectoring_longitude_angle_deg, thrustvectoring_latitude_angle_deg)


        template<typename StatesType, typename T>
        static constexpr states_type<T> derivatives(const StatesType &x,
                                                    const std::array<T, _num_states> &actuator_demands_deg,
                                                    const std::array<T, _num_states> &taus_s = { T{ default_parameters::tau_stabilator_s },
                                                                                                 T{ default_parameters::tau_leading_edge_flap_s },
                                                                                                 T{ default_parameters::tau_speedbreak_s },
                                                                                                 T{ default_parameters::tau_aileron_s },
                                                                                                 T{ default_parameters::tau_rudder_s },
                                                                                                 T{ default_parameters::tau_thrustvectoring_longitude_angle_s },
                                                                                                 T{ default_parameters::tau_thrustvectoring_latitude_angle_s } },
                                                    const std::array<T, _num_states> &maxabs_rates_degps = { T{ default_parameters::maxabs_stabilator_rate_degps },
                                                                                                             T{ default_parameters::maxabs_leading_edge_flap_rate_degps },
                                                                                                             T{ default_parameters::maxabs_speedbreak_rate_degps },
                                                                                                             T{ default_parameters::maxabs_aileron_rate_degps },
                                                                                                             T{ default_parameters::maxabs_rudder_rate_degps },
                                                                                                             T{ default_parameters::maxabs_thrustvectoring_longitude_angle_rate_degps },
                                                                                                             T{ default_parameters::maxabs_thrustvectoring_latitude_angle_rate_degps } });

        template<typename StatesType,
                 typename T = tr7::typeof_bracket_operator<StatesType> >
        static constexpr states_type<T> state_limiters(const StatesType &x,
                                                       const std::array<T, _num_states> &max_deflections_deg = { T{ default_parameters::max_stabilator_deg },
                                                                                                                 T{ default_parameters::max_leading_edge_flap_deg },
                                                                                                                 T{ default_parameters::max_speedbreak_deg },
                                                                                                                 T{ default_parameters::max_aileron_deg },
                                                                                                                 T{ default_parameters::max_rudder_deg },
                                                                                                                 T{ default_parameters::max_thrustvectoring_longitude_angle_deg },
                                                                                                                 T{ default_parameters::max_thrustvectoring_latitude_angle_deg } },
                                                       const std::array<T, _num_states> &min_deflections_deg = { T{ default_parameters::min_stabilator_deg },
                                                                                                                 T{ default_parameters::min_leading_edge_flap_deg },
                                                                                                                 T{ default_parameters::min_speedbreak_deg },
                                                                                                                 T{ default_parameters::min_aileron_deg },
                                                                                                                 T{ default_parameters::min_rudder_deg },
                                                                                                                 T{ default_parameters::min_thrustvectoring_longitude_angle_deg },
                                                                                                                 T{ default_parameters::min_thrustvectoring_latitude_angle_deg } });

    };


    //
    // "F16_Nguyen_actuators" impl.
    //

    template<typename StatesType, typename T>
    constexpr auto F16_Nguyen_actuators::derivatives(const StatesType &x,
                                                     const std::array<T, _num_states> &actuator_demands_deg,
                                                     const std::array<T, _num_states> &taus_s,
                                                     const std::array<T, _num_states> &maxabs_rates_degps) -> states_type<T>
    {
        //
        // Evaluates the model's state equation, i.e.,
        // "xdot = F(x, y, u, ...)", and returns the
        // model's "xdot".
        //

        ensure_minimum_size<_num_states>(x);
        states_type<T> xdot{};


        // calculate a state equation per actuator,
        // they're all the same
        tr7::tuple_index_for(xdot, [&](auto &xdot_i, auto i)
                                   {
                                       const auto ddot_degps = (actuator_demands_deg[i] - x[i]) / taus_s[i];

                                       if (tr7::AD_abs(ddot_degps) <= maxabs_rates_degps[i]) {
                                           xdot_i = ddot_degps;
                                       }
                                       else if (ddot_degps >= T{ 0 }) {
                                           xdot_i = maxabs_rates_degps[i];
                                       }
                                       else {
                                           xdot_i = -maxabs_rates_degps[i];
                                       }

                                   });


        // return the model's statedots
        return xdot;

    }


    template<typename StatesType,
             typename T>
    constexpr auto F16_Nguyen_actuators::state_limiters(const StatesType &x,
                                                        const std::array<T, _num_states> &max_deflections_deg,
                                                        const std::array<T, _num_states> &min_deflections_deg) -> states_type<T>
    {
        //
        // Evaluates the model's state limiter equation, i.e.,
        // "x = SL(x, ...)", and returns the model's limited "x".
        //

        ensure_minimum_size<_num_states>(x);
        states_type<T> xlim{};


        // limit each actuator deflection to its admissible range
        tr7::tuple_index_for(xlim, [&](auto &xlim_i, auto i)
                                   {
                                       if (x[i] < min_deflections_deg[i]) {
                                           xlim_i = min_deflections_deg[i];
                                       }
                                       else if (x[i] > max_deflections_deg[i]) {
                                           xlim_i = max_deflections_deg[i];
                                       }
                                       else {
                                           xlim_i = x[i];
                                       }

                                   });


        // return the limited states
        return xlim;

    }


}


#endif