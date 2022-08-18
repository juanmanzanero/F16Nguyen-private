#ifndef DYNAMIC_MODEL_F16_H
#define DYNAMIC_MODEL_F16_H

#include "F16_Nguyen/F16_Nguyen_plant.h"
#include "lion/thirdparty/include/cppad/cppad.hpp"
#include "lion/io/Xml_document.h"
#include "lion/foundation/types.h"
#include "src/core/vehicles/road_curvilinear.h"
#include "lion/math/matrix_extensions.h"
#include "F16_Nguyen/F16_Nguyen_trim.h"

#include "src/core/applications/optimal_laptime.h"
#include "/home/r660/r660391/storage/fastest-lap/src/core/vehicles/dynamic_model.h"


template<typename Track_t = Track_by_polynomial>
class F16_optimal_control : public Dynamic_model<CppAD::AD<double>>
{
 public:
    constexpr static size_t IX = F16_Nguyen::F16_Nguyen_plant::state_names::xEarth_m;
    constexpr static size_t IY = F16_Nguyen::F16_Nguyen_plant::state_names::yEarth_m;

    // Just to make sure...
    static_assert(IX == 10);
    static_assert(IY == 11);

    //! These are required by optimal_control.h
    constexpr static size_t ITIME = IX;
    constexpr static size_t IN = IY;

    using Track_type = Track_t;
    using Timeseries_type = CppAD::AD<double>;
    using Road_type = Road_curvilinear<Timeseries_type,Track_t, ITIME, 0>;
    using Dynamic_model_t = F16_optimal_control;

    //! TODO this is temporal, it is only required to compute the sensitivity analysis to laptime
    struct Chassis_type
    {
        constexpr static size_t IU = 0;
        constexpr static size_t IV = 0;
    };

    F16_optimal_control(const std::string& aeropath, const std::string& enginepath, Track_t& track, scalar mass, scalar xcg_per_MAC) : _pl(aeropath, enginepath), _road(track), _mass(mass), _xcg_per_MAC(xcg_per_MAC) {}


    //! The number of state variables
    constexpr static size_t NSTATE    = F16_Nguyen::F16_Nguyen_plant::num_states() ;

    //! The number of algebraic variables
    constexpr static size_t NALGEBRAIC = 0;

    //! The number of control variables
    constexpr static size_t NCONTROL  = 8 ; 

    //! Operator() 
    std::pair<std::array<Timeseries_type,NSTATE>,std::array<Timeseries_type,NALGEBRAIC>> operator()(const std::array<Timeseries_type,NSTATE>& q,
                                                                                                          const std::array<Timeseries_type,NALGEBRAIC>& qa,
                                                                                                          const std::array<Timeseries_type,NCONTROL>& u,
                                                                                                          scalar t)
    {
        // (1) Get the plant default inputs: get them scalar, convert to CppAD
        auto plant_inputs_scalar = _pl.default_inputs();
        plant_inputs_scalar[F16_Nguyen::F16_Nguyen_plant::input_names::mass_kg]     = _mass;
        plant_inputs_scalar[F16_Nguyen::F16_Nguyen_plant::input_names::xcg_per_MAC] = _xcg_per_MAC;

        std::array<Timeseries_type,std::tuple_size<decltype(plant_inputs_scalar)>::value> plant_inputs;
        std::copy(plant_inputs_scalar.cbegin(), plant_inputs_scalar.cend(), plant_inputs.begin());

        // (2) Update the inputs with the control variables: (only the NCONTROL first ones)
        std::copy(u.cbegin(), u.cend(), plant_inputs.begin());

        // (3) Compute plant outputs: assumption: the plant does not depend on (x,y). They will not be (x,y) but (time,n)
        const auto y = _pl.outputs(q,plant_inputs);

        // (4) Compute plant time derivative: the road will use psi=0, and the 'earth' velocities
        auto dqdt = _pl.derivatives(q,y,plant_inputs);
    
        // (5) Compute road time derivative -------------------------------------------------:-

        // (5.1) Update road: computes Frenet frame at the given arclength
        u_earth = dqdt[IX];
        v_earth = dqdt[IY];
        std::array<Timeseries_type,ITIME+3> q_road{0.0};

        theta = track_heading_angle_at(t);
        q_road[ITIME] = q[ITIME];
        q_road[IN]    = q[IN];
        q_road[IN+1]  = -theta;

        _road.set_state_and_controls(t, q_road, std::array<Timeseries_type,0>{});

        // (5.2) Call update for road: this computes road variables time derivative. We will bypass the rotations, which are left for the aircraft
        _road.update(u_earth, v_earth, 0.0);

        // (5.3) Download the road time derivative vector
        std::array<Timeseries_type,ITIME+3> dqdt_road;
        _road.get_state_derivative(dqdt_road);

        // (5.4) Fill the results into the global dqdt
        dqdt[IX] = dqdt_road[ITIME]; dqdt[IY] = dqdt_road[IN];

        // (6) Apply the chain rule to the entire dqdt: dqds = dqdt.dtds
        for (auto it = dqdt.begin(); it != dqdt.end(); ++it)
            (*it) *= _road.get_dtimedt();

        // (7) Save some outputs for later use
        aoa = y[F16_Nguyen::F16_Nguyen_plant::output_names::aoa_deg];
        aos = y[F16_Nguyen::F16_Nguyen_plant::output_names::aos_deg];
        yaw_deg = y[F16_Nguyen::F16_Nguyen_plant::output_names::yaw_deg];
        roll_deg = y[F16_Nguyen::F16_Nguyen_plant::output_names::roll_deg];
        gamma_deg = y[F16_Nguyen::F16_Nguyen_plant::output_names::ground_flight_path_angle_deg];
        nz_g = y[F16_Nguyen::F16_Nguyen_plant::output_names::Nz_g];

    //    std::cout << yaw_deg << " " << y[F16_Nguyen::F16_Nguyen_plant::output_names::pitch_deg] << " " << gamma_deg << std::endl;

        return {dqdt, {}};
    }

    const scalar track_heading_angle_at(const scalar t) 
    {
        _road.update_track(t);
        return  _road.get_heading_angle();
    }

    F16_Nguyen::F16_Nguyen_plant& get_plant() { return _pl; }

    Road_type& get_road() { return _road; }

    const Road_type& get_road() const { return _road; }

    // Optimal control extra definitions -------------------------------------------:-

    static constexpr const size_t N_OL_EXTRA_CONSTRAINTS = 6;

    std::tuple<std::vector<scalar>,std::vector<scalar>> optimal_laptime_derivative_control_bounds() const
    {   
        throw std::runtime_error("[ERROR] std::tuple<std::vector<scalar>,std::vector<scalar>> optimal_laptime_derivative_control_bounds() const ->"
                                 " not implemented");
        return {{}, {}};
    }

    std::pair<std::vector<scalar>,std::vector<scalar>> optimal_laptime_extra_constraints_bounds(const scalar s) const
    {   
        return
        {   
//          {-15.0,-15.0},
//          { 15.0, 15.0}
            {-20.0,-20.0,0.5,-90.0, -10.0, -4.0},
            { 20.0, 20.0,2.0, 90.0,  10.0,  10.0}
        };
    }

    std::array<Timeseries_type,N_OL_EXTRA_CONSTRAINTS> optimal_laptime_extra_constraints() const
    {
        const auto cos_angle = cos(yaw_deg*DEG - theta);

//      return {aoa, aos};
        return {aoa, aos, cos_angle, roll_deg, gamma_deg, nz_g};
    }

     
    //! Get state and control upper and lower values
    struct State_and_control_upper_lower_and_default_values
    {
        std::array<scalar,NSTATE> q_def;
        std::array<scalar,NSTATE> q_lb;
        std::array<scalar,NSTATE> q_ub;
        std::array<scalar,NALGEBRAIC> qa_def;
        std::array<scalar,NALGEBRAIC> qa_lb;
        std::array<scalar,NALGEBRAIC> qa_ub;
        std::array<scalar,NCONTROL> u_def;
        std::array<scalar,NCONTROL> u_lb;
        std::array<scalar,NCONTROL> u_ub;
    };

    static std::tuple<std::string,std::array<std::string,NSTATE>,std::array<std::string,NALGEBRAIC>,std::array<std::string,NCONTROL>>
        get_state_and_control_names() 
    {
        return 
        { 
            "arclength",
            { "u_mps","v_mps","w_mps","p_radps","q_radps","r_radps","qw_body2Earth","qx_body2Earth","qy_body2Earth","qz_body2Earth","xEarth_m","yEarth_m","zEarth_m","dh_deg","dlef_deg","dsb_deg","da_deg","dr_deg","thrustvectoring_longitude_angle_deg","thrustvectoring_latitude_angle_deg","P3_percent"},
            {},
            {"dh_dmd_deg","dlef_dmd_deg","dsb_dmd_deg","da_dmd_deg","dr_dmd_deg","thrustvectoring_longitude_angle_dmd_deg","thrustvectoring_latitude_angle_dmd_deg","throttle_percent"}
        };
    }
  
    State_and_control_upper_lower_and_default_values get_state_and_control_upper_lower_and_default_values() const
    {
        // (1) Define outputs
        std::array<scalar,NSTATE> q_def;        
        std::array<scalar,NSTATE> q_lb;
        std::array<scalar,NSTATE> q_ub;
        std::array<scalar,NALGEBRAIC> qa_def;
        std::array<scalar,NALGEBRAIC> qa_lb;
        std::array<scalar,NALGEBRAIC> qa_ub;
        std::array<scalar,NCONTROL> u_def;
        std::array<scalar,NCONTROL> u_lb;
        std::array<scalar,NCONTROL> u_ub;

        // (2) State
        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::u_mps] = 100.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::u_mps] = 1000.0;
        
        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::v_mps] = -100.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::v_mps] =  100.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::w_mps] = -100.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::w_mps] =  100.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::p_radps] = -10.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::p_radps] =  10.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::q_radps] = -10.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::q_radps] =  10.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::r_radps] = -10.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::r_radps] =  10.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qw_body2Earth] = -1.001;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qw_body2Earth] = 1.001;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qx_body2Earth] = -1.001;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qx_body2Earth] = 1.001;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qy_body2Earth] = -1.001;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qy_body2Earth] = 1.001;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qz_body2Earth] = -1.001;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qz_body2Earth] = 1.001;

        // (This is time, not x)
        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::xEarth_m] = 0.0; 
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::xEarth_m] = std::numeric_limits<double>::max();        

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::yEarth_m] = std::numeric_limits<double>::lowest();
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::yEarth_m] = std::numeric_limits<double>::max();        

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::zEarth_m] = -100.0; 
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::zEarth_m] =  10.0;    

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dh_deg] = -10.0; 
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dh_deg] =  10.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dlef_deg] =  -5.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dlef_deg] =   5.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dsb_deg] =  -5.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dsb_deg] =   5.0;
    
        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::da_deg] = -15.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::da_deg] =  15.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dr_deg] = -15.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dr_deg] =  15.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_longitude_angle_deg] = -5.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_longitude_angle_deg] =  5.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_latitude_angle_deg] = -5.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_latitude_angle_deg] =  5.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::P3_percent] = 0.00;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::P3_percent] = 80.0;

        // (3) Control vector
        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg] = dh_dmd_deg_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg] = -10.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg] = 10.0;

        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg] = dlef_dmd_deg_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg] = 15.0;

        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::dsb_dmd_deg] = dsb_dmd_deg_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dsb_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dsb_dmd_deg] = 15.0;
    
        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::da_dmd_deg] = da_dmd_deg_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::da_dmd_deg] = -10.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::da_dmd_deg] = 10.0;

        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::dr_dmd_deg] = dr_dmd_deg_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dr_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dr_dmd_deg] = 15.0;

        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_longitude_angle_dmd_deg] = thrustvectoring_longitude_angle_dmd_deg_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_longitude_angle_dmd_deg] = -5.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_longitude_angle_dmd_deg] =  5.0;

        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_latitude_angle_dmd_deg] = thrustvectoring_latitude_angle_dmd_deg_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_latitude_angle_dmd_deg] = -5.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_latitude_angle_dmd_deg] =  5.0;

        u_def[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent] = throttle_percent_default;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent] = 0.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent] = 80.0;

        // (4) Return
        return (State_and_control_upper_lower_and_default_values)
        {
            .q_def  = q_def , .q_lb  = q_lb , .q_ub  = q_ub ,
            .qa_def = qa_def, .qa_lb = qa_lb, .qa_ub = qa_ub,
            .u_def  = u_def , .u_lb  = u_lb , .u_ub  = u_ub
        };
    }
  
    Timeseries_type aoa;
    Timeseries_type aos;
    Timeseries_type theta;
    Timeseries_type u_earth;
    Timeseries_type v_earth;
    Timeseries_type yaw_deg;
    Timeseries_type roll_deg;
    Timeseries_type gamma_deg;
    Timeseries_type nz_g;

    scalar dh_dmd_deg_default                              = 0.0;
    scalar dlef_dmd_deg_default                            = 0.0;
    scalar dsb_dmd_deg_default                             = 0.0;
    scalar da_dmd_deg_default                              = 0.0;
    scalar dr_dmd_deg_default                              = 0.0;
    scalar thrustvectoring_longitude_angle_dmd_deg_default = 0.0;
    scalar thrustvectoring_latitude_angle_dmd_deg_default  = 0.0;
    scalar throttle_percent_default                        = 10.0;

 private:
    F16_Nguyen::F16_Nguyen_plant _pl;
    scalar _mass;
    scalar _xcg_per_MAC;

    Road_type _road;    //! 0,0 is preferred. I will just pass a size-3 array and manually place stuff in the big one

};

template<typename T>
std::pair<F16_Nguyen::F16_Nguyen_plant::states_type<T>,std::array<T,F16_optimal_control<>::NCONTROL>> trim2state
    (const F16_Nguyen::F16_Nguyen_trim::outputs_type<T>& trim_output, const F16_Nguyen::F16_Nguyen_trim::inputs_type<T>& trim_inputs)
{
    auto q = F16_Nguyen::F16_Nguyen_plant::states_type<T>{};
    auto u = std::array<T,F16_optimal_control<>::NCONTROL>{0.0};
    
    // Trim: pitch_deg, dh_deg, aoa_deg, dlef_deg, P3_percent, throttle_percent
    
    // Set pitch: quaternion is q = cos(pitch/2) + j.sin(pitch/2)
    const auto& pitch_deg = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::pitch_deg];

    q[F16_Nguyen::F16_Nguyen_plant::state_names::qw_body2Earth] = cos(0.5*pitch_deg*DEG);
    q[F16_Nguyen::F16_Nguyen_plant::state_names::qx_body2Earth] = 0.0;
    q[F16_Nguyen::F16_Nguyen_plant::state_names::qy_body2Earth] = sin(0.5*pitch_deg*DEG);
    q[F16_Nguyen::F16_Nguyen_plant::state_names::qz_body2Earth] = 0.0;

    // Set dh_deg and dh_dmd_deg
    q[F16_Nguyen::F16_Nguyen_plant::state_names::dh_deg] = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::dh_deg];
    u[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg] = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::dh_deg];

    // Set velocity
    const auto ZP_m = tr7::ft2m(trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::ZP_ft]);
    const auto TAS_mps = F16_Nguyen::cas2tas(tr7::kn2mps(trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::KCAS]),ZP_m);
    q[F16_Nguyen::F16_Nguyen_plant::state_names::u_mps] = TAS_mps*cos(trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::aoa_deg]*DEG);
    q[F16_Nguyen::F16_Nguyen_plant::state_names::w_mps] = TAS_mps*sin(trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::aoa_deg]*DEG);

    // Set altitude
    q[F16_Nguyen::F16_Nguyen_plant::state_names::zEarth_m] = -ZP_m;

    // Set leading edge flap
    q[F16_Nguyen::F16_Nguyen_plant::state_names::dlef_deg] = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::dlef_deg];
    u[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg] = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::dlef_deg];

    // Set P3 percent
    q[F16_Nguyen::F16_Nguyen_plant::state_names::P3_percent] = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::P3_percent];

    // Set throttle percent
    u[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent] = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::throttle_percent];

    return {q, u};
}

#endif
