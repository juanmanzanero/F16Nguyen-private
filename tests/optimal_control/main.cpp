#include <iomanip>

#include "F16_Nguyen/F16_Nguyen_plant.h"
#include "lion/thirdparty/include/cppad/cppad.hpp"
#include "lion/io/Xml_document.h"
#include "lion/foundation/types.h"
#include "/Users/juanmanzanero/Documents/software/fastest-lap/src/core/vehicles/road_curvilinear.h"
#include "/Users/juanmanzanero/Documents/software/fastest-lap/src/core/vehicles/track_by_arcs.h"
#include "lion/math/matrix_extensions.h"
#include "F16_Nguyen/F16_Nguyen_trim.h"

#include "/Users/juanmanzanero/Documents/software/fastest-lap/src/core/applications/optimal_laptime.h"


class F16_optimal_control
{
 public:
    using Timeseries_type = CppAD::AD<double>;
    using Road_type = Road_curvilinear<Timeseries_type,Track_by_arcs, 0, 0>;
    using Dynamic_model_t = F16_optimal_control;

    F16_optimal_control(const std::string& aeropath, const std::string& enginepath, Track_by_arcs& track, scalar mass, scalar xcg_per_MAC) : _pl(aeropath, enginepath), _road(track), _mass(mass), _xcg_per_MAC(xcg_per_MAC) {}

    constexpr static size_t IX = F16_Nguyen::F16_Nguyen_plant::state_names::xEarth_m;
    constexpr static size_t IY = F16_Nguyen::F16_Nguyen_plant::state_names::yEarth_m;

    // Just to make sure...
    static_assert(IX == 10);
    static_assert(IY == 11);

    //! These are required by optimal_control.h
    constexpr static size_t ITIME = IX;
    constexpr static size_t IN = IY;


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
        const auto u_earth = dqdt[IX];
        const auto v_earth = dqdt[IY];
        const auto psi     = 0.0;
        std::array<Timeseries_type,3> q_road = {q[ITIME], q[IN], 0.0};

        _road.set_state_and_controls(t, q_road, std::array<Timeseries_type,0>{});

        // (5.2) Call update for road: this computes road variables time derivative. We will bypass the rotations, which are left for the aircraft
        _road.update(u_earth, v_earth, 0.0);

        // (5.3) Download the road time derivative vector
        std::array<Timeseries_type,3> dqdt_road;
        _road.get_state_derivative(dqdt_road);

        // (5.4) Fill the results into the global dqdt
        dqdt[IX] = dqdt_road[ITIME]; dqdt[IY] = dqdt_road[IN];

        // (6) Apply the chain rule to the entire dqdt: dqds = dqdt.dtds
        for (auto it = dqdt.begin(); it != dqdt.end(); ++it)
            (*it) *= _road.get_dtimedt();

        return {dqdt, {}};
    }

    F16_Nguyen::F16_Nguyen_plant& get_plant() { return _pl; }

    Road_type& get_road() { return _road; }

    const Road_type& get_road() const { return _road; }

    // Optimal control extra definitions -------------------------------------------:-

    static constexpr const size_t N_OL_EXTRA_CONSTRAINTS = 0;

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
            {},
            {}
        };
    }

    std::array<Timeseries_type,N_OL_EXTRA_CONSTRAINTS> optimal_laptime_extra_constraints() const
    {
        return {};
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

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qw_body2Earth] = 0.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qw_body2Earth] = 1.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qx_body2Earth] = 0.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qx_body2Earth] = 1.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qy_body2Earth] = 0.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qy_body2Earth] = 1.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::qz_body2Earth] = 0.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::qz_body2Earth] = 1.0;

        // (This is time, not x)
        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::xEarth_m] = 0.0; 
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::xEarth_m] = std::numeric_limits<double>::max();        

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::yEarth_m] = std::numeric_limits<double>::lowest();
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::yEarth_m] = std::numeric_limits<double>::max();        

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::zEarth_m] = -100.0; 
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::zEarth_m] =  1.0;    

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dh_deg] = -15.0; 
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dh_deg] =  15.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dlef_deg] =  -15.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dlef_deg] =   15.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dsb_deg] =  -15.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dsb_deg] =   15.0;
    
        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::da_deg] = -15.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::da_deg] =  15.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::dr_deg] = -15.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::dr_deg] =  15.0;


        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_longitude_angle_deg] = -1.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_longitude_angle_deg] =  1.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_latitude_angle_deg] = -1.0;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::thrustvectoring_latitude_angle_deg] =  1.0;

        q_lb[F16_Nguyen::F16_Nguyen_plant::state_names::P3_percent] = 0.01;
        q_ub[F16_Nguyen::F16_Nguyen_plant::state_names::P3_percent] = 100.0;

        // (3) Control vector
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg] = 15.0;

        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg] = 15.0;

        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dsb_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dsb_dmd_deg] = 15.0;
    
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::da_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::da_dmd_deg] = 15.0;

        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::dr_dmd_deg] = -15.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::dr_dmd_deg] = 15.0;

        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_longitude_angle_dmd_deg] = -5.0;
        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_longitude_angle_dmd_deg] =  5.0;

        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_latitude_angle_dmd_deg] = -5.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_latitude_angle_dmd_deg] =  5.0;

        u_lb[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent] = 0.0;
        u_ub[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent] = 100.0;

        // (4) Return
        return (State_and_control_upper_lower_and_default_values)
        {
            .q_def  = q_def , .q_lb  = q_lb , .q_ub  = q_ub ,
            .qa_def = qa_def, .qa_lb = qa_lb, .qa_ub = qa_ub,
            .u_def  = u_def , .u_lb  = u_lb , .u_ub  = u_ub
        };
    }
  

 private:
    F16_Nguyen::F16_Nguyen_plant _pl;
    scalar _mass;
    scalar _xcg_per_MAC;

    Road_type _road;    //! 0,0 is preferred. I will just pass a size-3 array and manually place stuff in the big one

};

template<typename T>
std::pair<F16_Nguyen::F16_Nguyen_plant::states_type<T>,std::array<T,F16_optimal_control::NCONTROL>> trim2state
    (const F16_Nguyen::F16_Nguyen_trim::outputs_type<T>& trim_output, const F16_Nguyen::F16_Nguyen_trim::inputs_type<T>& trim_inputs)
{
    auto q = F16_Nguyen::F16_Nguyen_plant::states_type<T>{};
    auto u = std::array<T,F16_optimal_control::NCONTROL>{0.0};
    
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

// First naive test to Optimal Control of a F16
int main()
{
    // (1) Database paths
    const std::string aeropath = "../../datasets/aero";
    const std::string enginepath = "../../datasets/engine";
    const scalar mass = 9.0e3;
    const scalar xcg_per_MAC = 0.55;

    // (2) Construct track
    Xml_document straight_xml("./straight.xml", true);
    Track_by_arcs straight(straight_xml,false);

    const size_t n = 100;
    const auto s = linspace(0.0,straight.get_total_length(),n+1);

    // (3) Construct F16
    F16_optimal_control vehicle(aeropath, enginepath, straight, mass, xcg_per_MAC);

    // (4) Trim F16
    auto trim_inputs = F16_Nguyen::F16_Nguyen_trim::default_inputs();
    trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::KCAS] = 250.0;
    trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::mass_kg] = mass;
    trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::xcg_per_MAC] = xcg_per_MAC;
    
    auto trim_inflags = F16_Nguyen::F16_Nguyen_trim::default_inflags();
    trim_inflags[F16_Nguyen::F16_Nguyen_trim::inflag_names::disable_lef] = false;

    auto [trim_success, trim_output] = F16_Nguyen::F16_Nguyen_trim::trim_plant(vehicle.get_plant(), trim_inflags, trim_inputs); 

    auto [q,u] = trim2state(trim_output, trim_inputs);

    if ( !trim_success ) throw std::runtime_error("Trim did not succeed");

    // (5) Check trim point: transform to cppad + run operator()
    std::array<CppAD::AD<scalar>, F16_Nguyen::F16_Nguyen_plant::num_states()> q_cppad;
    std::array<CppAD::AD<scalar>, F16_optimal_control::NCONTROL> u_cppad;

    std::copy(q.cbegin(), q.cend(), q_cppad.begin());
    std::copy(u.cbegin(), u.cend(), u_cppad.begin());

    for (auto& dqdt_i : vehicle(q_cppad,{},u_cppad,0.0).first )
        assert(std::abs(dqdt_i) < 1.0e-13 );

    std::cout << "[info] Trim successful:" << std::endl;
    std::cout << "[info]     -> state: " << q << std::endl;
    std::cout << "[info]     -> controls: " << u << std::endl;

  
    // (6) Run optimal control in a straight line
    auto control_variables = typename Optimal_laptime<F16_optimal_control>::template Control_variables<>{};

    for (size_t i = 0; i < F16_optimal_control::NCONTROL; ++i)
        control_variables[i] = Optimal_laptime<F16_optimal_control>::create_full_mesh({n+1,u[i]}, 1.0e-2);
    
    Optimal_laptime<F16_optimal_control>::Options opts;
    opts.print_level = 5;
    Optimal_laptime opt_laptime(s, true, false, vehicle, {n+1,q}, {}, control_variables, opts);
   

}
