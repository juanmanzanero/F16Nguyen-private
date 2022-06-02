#ifndef F16_NGUYEN_F16_NGUYEN_PLANT_H
#define F16_NGUYEN_F16_NGUYEN_PLANT_H
#pragma once


#include "F16_Nguyen/F16_Nguyen_air_kinematics.h"
#include "F16_Nguyen/F16_Nguyen_aerodynamics.h"
#include "F16_Nguyen/F16_Nguyen_engine.h"
#include "F16_Nguyen/F16_Nguyen_flat_Earth_body_dynamics.h"
#include "F16_Nguyen/F16_Nguyen_air_dynamics.h"
#include "F16_Nguyen/F16_Nguyen_actuators.h"


//
// Defines the "F16_Nguyen_plant" class, that collects
// all of the submodels we've written into a single
// state-space, following reference "Simulator study of
// stall/post-stall characteristics of a fighter airplane
// with relaxed longitudinal static stability", L.T.Nguyen
// et al., NASA-TP-1538, 1979.
//


namespace F16_Nguyen
{
    class F16_Nguyen_plant
    {
        //
        // Represents the whole F16 model, a
        // concatenation of each submodel we have.
        //

    public:

        using data_type = double;


        // we'll provide an interface of inputs for this
        // "aggregate" model
        DECLARE_MODEL_INPUTS(dh_dmd_deg, dlef_dmd_deg, dsb_dmd_deg,
                             da_dmd_deg, dr_dmd_deg,
                             thrustvectoring_longitude_angle_dmd_deg, thrustvectoring_latitude_angle_dmd_deg,
                             throttle_percent,
                             aoadot_degps,
                             mass_kg,
                             Ixx_kgm2, Ixy_kgm2, Ixz_kgm2, Iyy_kgm2, Iyz_kgm2, Izz_kgm2,
                             thrust_point_of_applicationx_bodyaxes_m, thrust_point_of_applicationy_bodyaxes_m, thrust_point_of_applicationz_bodyaxes_m,
                             xcg_per_MAC, ycg_per_semiwingspan, MAC_m, wingspan_m, wing_surface_m2,
                             DISA_K)


        // now we declare the model's states & outputs, by concatenating the
        // states & outputs of the submodels
        using state_names = detail::glue_state_names<F16_Nguyen_flat_Earth_body_dynamics,
                                                     F16_Nguyen_actuators,
                                                     F16_Nguyen_engine>;

        static constexpr std::size_t _num_states = state_names::_num_names;
        static constexpr std::size_t num_states() { return _num_states; }

        template<typename T>
        using states_type = std::array<T, _num_states>;


        using output_names = detail::glue_output_names<F16_Nguyen_air_kinematics,
                                                       F16_Nguyen_engine,
                                                       F16_Nguyen_aerodynamics,
                                                       F16_Nguyen_flat_Earth_body_dynamics,
                                                       F16_Nguyen_air_dynamics>;

        static constexpr std::size_t _num_outputs = output_names::_num_names;
        static constexpr std::size_t num_outputs() { return _num_outputs; }

        template<typename T>
        using outputs_type = std::array<T, _num_outputs>;


        // helper aliases for the states: use the first one to access their [_begin, _end) ranges,
        // and the second one to access state names that are ambiguous amongst different models
        // (unambiguous names can be accessed simply via "state_names::")
        using states_body_dynamics = state_names::_root;
        using state_names_body_dynamics = states_body_dynamics::_names;
        static_assert(std::is_same_v<F16_Nguyen_flat_Earth_body_dynamics, states_body_dynamics::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model states!");

        using states_actuators = states_body_dynamics::_next;
        using state_names_actuators = states_actuators::_names;
        static_assert(std::is_same_v<F16_Nguyen_actuators, states_actuators::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model states!");

        using states_engine = states_actuators::_next;
        using state_names_engine = states_engine::_names;
        static_assert(std::is_same_v<F16_Nguyen_engine, states_engine::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model states!");


        // helper aliases for the outputs: use the first one to access their [_begin, _end) ranges,
        // and the second one to access output names that are ambiguous amongst different models
        // (unambiguous names can be accessed simply via "output_names::")
        using outputs_air_kinematics = output_names::_root;
        using output_names_air_kinematics = outputs_air_kinematics::_names;
        static_assert(std::is_same_v<F16_Nguyen_air_kinematics, outputs_air_kinematics::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model outputs!");

        using outputs_engine = outputs_air_kinematics::_next;
        using output_names_engine = outputs_engine::_names;
        static_assert(std::is_same_v<F16_Nguyen_engine, outputs_engine::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model outputs!");

        using outputs_aerodynamics = outputs_engine::_next;
        using output_names_aerodynamics = outputs_aerodynamics::_names;
        static_assert(std::is_same_v<F16_Nguyen_aerodynamics, outputs_aerodynamics::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model outputs!");

        using outputs_body_dynamics = outputs_aerodynamics::_next;
        using output_names_body_dynamics = outputs_body_dynamics::_names;
        static_assert(std::is_same_v<F16_Nguyen_flat_Earth_body_dynamics, outputs_body_dynamics::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model outputs!");

        using outputs_air_dynamics = outputs_body_dynamics::_next;
        using output_names_air_dynamics = outputs_air_dynamics::_names;
        static_assert(std::is_same_v<F16_Nguyen_air_dynamics, outputs_air_dynamics::_model>, "F16_Nguyen::F16_Nguyen_plant: inconsistent model outputs!");


        // export states & output names as arrays of std::strings
        static inline const std::array<std::string, _num_states>& state_names2str()
        {
            static const auto arr = detail::strarray_cat(state_names::_root::_model::state_names2str(),
                                                         state_names::_root::_next::_model::state_names2str(),
                                                         state_names::_root::_next::_next::_model::state_names2str());
            return arr;

        }

        static inline const std::array<std::string, _num_outputs>& output_names2str()
        {
            static const auto arr = detail::strarray_cat(output_names::_root::_model::output_names2str(),
                                                         output_names::_root::_next::_model::output_names2str(),
                                                         output_names::_root::_next::_next::_model::output_names2str(),
                                                         output_names::_root::_next::_next::_next::_model::output_names2str(),
                                                         output_names::_root::_next::_next::_next::_next::_model::output_names2str());
            return arr;

        }


        F16_Nguyen_plant() = default;

        template<typename... Args>
        F16_Nguyen_plant(Args&&... args)
        {
            build(std::forward<Args>(args)...);
        }

        void build(const std::string &path_aerodataset,
                   const std::string &path_enginedataset)
        {
            _aerodynamics.build(path_aerodataset);
            _engine.build(path_enginedataset);

        }


        const F16_Nguyen_aerodynamics& aerodynamics() const { return _aerodynamics; }
        const F16_Nguyen_engine& engine() const { return _engine; }


        template<typename T = data_type>
        static constexpr inputs_type<T> default_inputs();


        template<typename StatesType,
                 typename T = tr7::typeof_bracket_operator<StatesType> >
        constexpr outputs_type<T> outputs(const StatesType &x,
                                          const T &throttle_percent,
                                          const T &aoadot_degps,
                                          const T &mass_kg,
                                          const T &Ixx_kgm2, const T &Ixy_kgm2, const T &Ixz_kgm2, const T &Iyy_kgm2, const T &Iyz_kgm2, const T &Izz_kgm2,
                                          const T &thrust_point_of_applicationx_bodyaxes_m, const T &thrust_point_of_applicationy_bodyaxes_m, const T &thrust_point_of_applicationz_bodyaxes_m,
                                          const T &xcg_per_MAC, const T &ycg_per_semiwingspan, const T &MAC_m, const T &wingspan_m, const T &wing_surface_m2,
                                          const T &DISA_K) const;

        template<typename StatesType, typename InputsType>
        constexpr auto outputs(const StatesType &x,
                               const InputsType &u) const
        {
            ensure_minimum_size<_num_inputs>(u);
            return outputs(x,
                           u[input_names::throttle_percent],
                           u[input_names::aoadot_degps],
                           u[input_names::mass_kg],
                           u[input_names::Ixx_kgm2], u[input_names::Ixy_kgm2], u[input_names::Ixz_kgm2], u[input_names::Iyy_kgm2], u[input_names::Iyz_kgm2], u[input_names::Izz_kgm2],
                           u[input_names::thrust_point_of_applicationx_bodyaxes_m], u[input_names::thrust_point_of_applicationy_bodyaxes_m], u[input_names::thrust_point_of_applicationz_bodyaxes_m],
                           u[input_names::xcg_per_MAC], u[input_names::ycg_per_semiwingspan], u[input_names::MAC_m], u[input_names::wingspan_m], u[input_names::wing_surface_m2],
                           u[input_names::DISA_K]);

        }


        template<typename StatesType, typename OutputsType,typename T>
        static constexpr states_type<T> derivatives(const StatesType &x,
                                                    const OutputsType &y,
                                                    const T &dh_dmd_deg, const T &dlef_dmd_deg, const T &dsb_dmd_deg,
                                                    const T &da_dmd_deg, const T &dr_dmd_deg,
                                                    const T &thrustvectoring_longitude_angle_dmd_deg, const T &thrustvectoring_latitude_angle_dmd_deg);

        template<typename StatesType, typename OutputsType, typename InputsType>
        static constexpr auto derivatives(const StatesType &x,
                                          const OutputsType &y,
                                          const InputsType &u)
        {
            ensure_minimum_size<_num_inputs>(u);
            return derivatives(x, y,
                               u[input_names::dh_dmd_deg], u[input_names::dlef_dmd_deg], u[input_names::dsb_dmd_deg],
                               u[input_names::da_dmd_deg], u[input_names::dr_dmd_deg],
                               u[input_names::thrustvectoring_longitude_angle_dmd_deg], u[input_names::thrustvectoring_latitude_angle_dmd_deg]);

        }


        template<typename StatesType,
                 typename T = tr7::typeof_bracket_operator<StatesType>,
                 typename... Args>
        static constexpr states_type<T> state_limiters(const StatesType &x, Args&&...);


    private:

        F16_Nguyen_aerodynamics _aerodynamics;
        F16_Nguyen_engine _engine;

    };


    //
    // "F16_Nguyen_plant" impl.
    //

    template<typename T>
    constexpr auto F16_Nguyen_plant::default_inputs() -> inputs_type<T>
    {
        //
        // Returns an array of default inputs so that the user
        // doesn't need to set all of them every time they
        // call the plant.
        //

        constexpr auto dh_dmd_deg = T{ 0 };
        constexpr auto dlef_dmd_deg = T{ 0 };
        constexpr auto dsb_dmd_deg = T{ 0 };
        constexpr auto da_dmd_deg = T{ 0 };
        constexpr auto dr_dmd_deg = T{ 0 };
        constexpr auto thrustvectoring_longitude_angle_dmd_deg = T{ 0 };
        constexpr auto thrustvectoring_latitude_angle_dmd_deg = T{ 0 };
        constexpr auto throttle_percent = T{ 0 };
        constexpr auto aoadot_degps = T{ 0 };

        constexpr auto mass_kg = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::mass_kg };
        constexpr auto Ixx_kgm2 = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::Ixx_kgm2 };
        constexpr auto Ixy_kgm2 = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::Ixy_kgm2 };
        constexpr auto Ixz_kgm2 = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::Ixz_kgm2 };
        constexpr auto Iyy_kgm2 = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::Iyy_kgm2 };
        constexpr auto Iyz_kgm2 = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::Iyz_kgm2 };
        constexpr auto Izz_kgm2 = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::Izz_kgm2 };

        constexpr auto thrust_point_of_applicationx_bodyaxes_m = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::thrust_point_of_applicationx_bodyaxes_m };
        constexpr auto thrust_point_of_applicationy_bodyaxes_m = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::thrust_point_of_applicationy_bodyaxes_m };
        constexpr auto thrust_point_of_applicationz_bodyaxes_m = T{ F16_Nguyen_flat_Earth_body_dynamics::default_parameters::thrust_point_of_applicationz_bodyaxes_m };

        constexpr auto xcg_per_MAC = T{ F16_Nguyen_aerodynamics::default_parameters::reference_xcg_per_MAC };
        constexpr auto ycg_per_semiwingspan = T{ F16_Nguyen_aerodynamics::default_parameters::reference_ycg_per_semiwingspan };
        constexpr auto MAC_m = T{ F16_Nguyen_aerodynamics::default_parameters::MAC_m };
        constexpr auto wingspan_m = T{ F16_Nguyen_aerodynamics::default_parameters::wingspan_m };
        constexpr auto wing_surface_m2 = T{ F16_Nguyen_aerodynamics::default_parameters::wing_surface_m2 };

        constexpr auto DISA_K = T{ 0 };

        return inputs_type<T>{ dh_dmd_deg, dlef_dmd_deg, dsb_dmd_deg,
                               da_dmd_deg, dr_dmd_deg,
                               thrustvectoring_longitude_angle_dmd_deg, thrustvectoring_latitude_angle_dmd_deg,
                               throttle_percent,
                               aoadot_degps,
                               mass_kg,
                               Ixx_kgm2, Ixy_kgm2, Ixz_kgm2, Iyy_kgm2, Iyz_kgm2, Izz_kgm2,
                               thrust_point_of_applicationx_bodyaxes_m, thrust_point_of_applicationy_bodyaxes_m, thrust_point_of_applicationz_bodyaxes_m,
                               xcg_per_MAC, ycg_per_semiwingspan, MAC_m, wingspan_m, wing_surface_m2,
                               DISA_K };

    }


    template<typename StatesType, typename T>
    constexpr auto F16_Nguyen_plant::outputs(const StatesType &x,
                                             const T &throttle_percent,
                                             const T &aoadot_degps,
                                             const T &mass_kg,
                                             const T &Ixx_kgm2, const T &Ixy_kgm2, const T &Ixz_kgm2, const T &Iyy_kgm2, const T &Iyz_kgm2, const T &Izz_kgm2,
                                             const T &thrust_point_of_applicationx_bodyaxes_m, const T &thrust_point_of_applicationy_bodyaxes_m, const T &thrust_point_of_applicationz_bodyaxes_m,
                                             const T &xcg_per_MAC, const T &ycg_per_semiwingspan, const T &MAC_m, const T &wingspan_m, const T &wing_surface_m2,
                                             const T &DISA_K) const -> outputs_type<T>
    {
        //
        // Evaluates the model's output equation, i.e.,
        // "y = G(x, u, ...)", and returns the model's "y".
        //

        ensure_minimum_size<_num_states>(x);
        outputs_type<T> y{};


        // call air kinematics
        std::copy_n(F16_Nguyen_air_kinematics::outputs(-x[state_names::zEarth_m],
                                                       x[state_names::u_mps], x[state_names::v_mps], x[state_names::w_mps],
                                                       x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps],
                                                       x[state_names::qw_body2Earth], x[state_names::qx_body2Earth], x[state_names::qy_body2Earth], x[state_names::qz_body2Earth],
                                                       DISA_K).cbegin(),
                    F16_Nguyen_air_kinematics::_num_outputs, &std::get<outputs_air_kinematics::_begin>(y));


        // call engine
        std::copy_n(_engine.outputs(make_const_static_contiguous_range<decltype(_engine)::_num_states>(&x[states_engine::_begin]),
                                    throttle_percent,
                                    y[output_names_air_kinematics::Mach],
                                    tr7::ft2m(y[output_names_air_kinematics::ZP_ft])).cbegin(),
                    _engine._num_outputs, &std::get<outputs_engine::_begin>(y));


        // call aerodynamics
        std::copy_n(_aerodynamics.outputs(y[output_names_air_kinematics::dynamic_pressure_Pa],
                                          y[output_names_air_kinematics::TAS_mps], y[output_names_air_kinematics::aoa_deg], y[output_names_air_kinematics::aos_deg],
                                          y[output_names_air_kinematics::paero_radps], y[output_names_air_kinematics::qaero_radps], y[output_names_air_kinematics::raero_radps],
                                          x[state_names::dh_deg], x[state_names::dlef_deg], x[state_names::dsb_deg],
                                          x[state_names::da_deg], x[state_names::dr_deg],
                                          tr7::deg2rad(aoadot_degps),
                                          xcg_per_MAC, ycg_per_semiwingspan, MAC_m, wingspan_m, wing_surface_m2).cbegin(),
                    _aerodynamics._num_outputs, &std::get<outputs_aerodynamics::_begin>(y));


        // call body dynamics
        std::copy_n(F16_Nguyen_flat_Earth_body_dynamics::outputs(make_const_static_contiguous_range<F16_Nguyen_flat_Earth_body_dynamics::_num_states>(&x[states_body_dynamics::_begin]),
                                                                 y[output_names_aerodynamics::Faerox_N], y[output_names_aerodynamics::Faeroy_N], y[output_names_aerodynamics::Faeroz_N],
                                                                 y[output_names_aerodynamics::Maerox_Nm], y[output_names_aerodynamics::Maeroy_Nm], y[output_names_aerodynamics::Maeroz_Nm],
                                                                 y[output_names_engine::engine_thrust_N], y[output_names_engine::engine_angular_momentum_kgm2ps],
                                                                 x[state_names::thrustvectoring_longitude_angle_deg], x[state_names::thrustvectoring_latitude_angle_deg],
                                                                 mass_kg,
                                                                 Ixx_kgm2, Ixy_kgm2, Ixz_kgm2, Iyy_kgm2, Iyz_kgm2, Izz_kgm2,
                                                                 thrust_point_of_applicationx_bodyaxes_m, thrust_point_of_applicationy_bodyaxes_m, thrust_point_of_applicationz_bodyaxes_m).cbegin(),
                    F16_Nguyen_flat_Earth_body_dynamics::_num_outputs, &std::get<outputs_body_dynamics::_begin>(y));


        // call air dynamics
        std::copy_n(F16_Nguyen_air_dynamics::outputs(y[output_names_body_dynamics::udot_mps2], y[output_names_body_dynamics::vdot_mps2], y[output_names_body_dynamics::wdot_mps2],
                                                     y[output_names_air_kinematics::VTASx_mps], y[output_names_air_kinematics::VTASy_mps], y[output_names_air_kinematics::VTASz_mps],
                                                     y[output_names_air_kinematics::TAS_mps], y[output_names_air_kinematics::aoa_deg], y[output_names_air_kinematics::aos_deg],
                                                     y[output_names_air_kinematics::true_heading_angle_deg], y[output_names_air_kinematics::flight_path_angle_deg], y[output_names_air_kinematics::flight_roll_angle_deg],
                                                     y[output_names_body_dynamics::Nx_g], y[output_names_body_dynamics::Ny_g], y[output_names_body_dynamics::Nz_g],
                                                     x[state_names::p_radps], x[state_names::q_radps], x[state_names::r_radps]).cbegin(),
                    F16_Nguyen_air_dynamics::_num_outputs, &std::get<outputs_air_dynamics::_begin>(y));


        // return the model's outputs
        return y;

    }


    template<typename StatesType, typename OutputsType, typename T>
    constexpr auto F16_Nguyen_plant::derivatives(const StatesType &x,
                                                 const OutputsType &y,
                                                 const T &dh_dmd_deg,
                                                 const T &dlef_dmd_deg,
                                                 const T &dsb_dmd_deg,
                                                 const T &da_dmd_deg,
                                                 const T &dr_dmd_deg,
                                                 const T &thrustvectoring_longitude_angle_dmd_deg,
                                                 const T &thrustvectoring_latitude_angle_dmd_deg) -> states_type<T>
    {
        //
        // Evaluates the model's state equation, i.e.,
        // "xdot = F(x, y, u, ...)", and returns the
        // model's "xdot".
        //

        ensure_minimum_size<_num_states>(x);
        ensure_minimum_size<_num_outputs>(y);
        states_type<T> xdot{};


        // call body dynamics
        std::copy_n(F16_Nguyen_flat_Earth_body_dynamics::derivatives(make_const_static_contiguous_range<F16_Nguyen_flat_Earth_body_dynamics::_num_states>(&x[states_body_dynamics::_begin]),
                                                                     make_const_static_contiguous_range<F16_Nguyen_flat_Earth_body_dynamics::_num_outputs>(&y[outputs_body_dynamics::_begin])).cbegin(),
                    F16_Nguyen_flat_Earth_body_dynamics::_num_states, &std::get<states_body_dynamics::_begin>(xdot));


        // call actuators
        std::copy_n(F16_Nguyen_actuators::derivatives(make_const_static_contiguous_range<F16_Nguyen_actuators::_num_states>(&x[states_actuators::_begin]),
                                                      std::array<T, states_actuators::_model::_num_states>{ dh_dmd_deg, dlef_dmd_deg, dsb_dmd_deg,
                                                                                                            da_dmd_deg, dr_dmd_deg,
                                                                                                            thrustvectoring_longitude_angle_dmd_deg,
                                                                                                            thrustvectoring_latitude_angle_dmd_deg }).cbegin(),
                    F16_Nguyen_actuators::_num_states, &std::get<states_actuators::_begin>(xdot));


        // call engine
        std::copy_n(F16_Nguyen_engine::derivatives(make_const_static_contiguous_range<F16_Nguyen_engine::_num_states>(&x[states_engine::_begin]),
                                                   make_const_static_contiguous_range<F16_Nguyen_engine::_num_outputs>(&y[outputs_engine::_begin])).cbegin(),
                    F16_Nguyen_engine::_num_states, &std::get<states_engine::_begin>(xdot));


        // return the model's statedots
        return xdot;

    }


    template<typename StatesType, typename T,
             typename... Args>
    constexpr auto F16_Nguyen_plant::state_limiters(const StatesType &x, Args&&...) -> states_type<T>
    {
        //
        // Evaluates the model's state limiter equation, i.e.,
        // "x = SL(x, ...)", and returns the model's limited "x".
        //

        ensure_minimum_size<_num_states>(x);
        states_type<T> xlim{};


        // call body dynamics
        std::copy_n(F16_Nguyen_flat_Earth_body_dynamics::state_limiters(
                    make_const_static_contiguous_range<F16_Nguyen_flat_Earth_body_dynamics::_num_states>(&x[states_body_dynamics::_begin])).cbegin(),
                    F16_Nguyen_flat_Earth_body_dynamics::_num_states, &std::get<states_body_dynamics::_begin>(xlim));


        // call actuators
        std::copy_n(F16_Nguyen_actuators::state_limiters(
                    make_const_static_contiguous_range<F16_Nguyen_actuators::_num_states>(&x[states_actuators::_begin])).cbegin(),
                    F16_Nguyen_actuators::_num_states, &std::get<states_actuators::_begin>(xlim));


        // call engine
        std::copy_n(F16_Nguyen_engine::state_limiters(
                    make_const_static_contiguous_range<F16_Nguyen_engine::_num_states>(&x[states_engine::_begin])).cbegin(),
                    F16_Nguyen_engine::_num_states, &std::get<states_engine::_begin>(xlim));


        // return the model's limited states
        return xlim;

    }


}


#endif