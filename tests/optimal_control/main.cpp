#include <iomanip>

#include "src/core/vehicles/track_by_arcs.h"
#include "src/core/vehicles/track_by_polynomial.h"
#include "dynamic_model_F16.h"

// First naive test to Optimal Control of a F16
int main()
{
    // (1) Database paths
    const std::string aeropath = "../../datasets/aero";
    const std::string enginepath = "../../datasets/engine";
    const scalar mass = 9.0e3;
    const scalar xcg_per_MAC = 0.55;

    // (2) Construct track
/*
    Xml_document catalunya_xml("/Users/juanmanzanero/Documents/software/fastest-lap/database/tracks/catalunya/catalunya_discrete.xml",true);
    Circuit_preprocessor catalunya_pproc(catalunya_xml);
    Track_by_polynomial catalunya(catalunya_pproc);
    const auto& s = catalunya_pproc.s;
    const auto& n = s.size();
*/

/*
    Xml_document ovaltrack_xml("./straight.xml",true);
    Track_by_arcs ovaltrack(ovaltrack_xml,1.0,false);
*/
  
    Xml_document ovaltrack_xml("./ovaltrack.xml",true);
    Track_by_arcs ovaltrack(ovaltrack_xml,60.0,true);
    
    constexpr const size_t n = 200;

    auto s = linspace(0.0,0.60*ovaltrack.get_total_length(), n+1);

    // (3) Construct F16
    F16_optimal_control vehicle(aeropath, enginepath, ovaltrack, mass, xcg_per_MAC);

    // (4) Trim F16
    auto trim_inputs = F16_Nguyen::F16_Nguyen_trim::default_inputs();
    trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::KCAS] = 400.0;
    trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::mass_kg] = mass;
    trim_inputs[F16_Nguyen::F16_Nguyen_trim::input_names::xcg_per_MAC] = xcg_per_MAC;
    
    auto trim_inflags = F16_Nguyen::F16_Nguyen_trim::default_inflags();
    trim_inflags[F16_Nguyen::F16_Nguyen_trim::inflag_names::disable_lef] = true;

    auto [trim_success, trim_output] = F16_Nguyen::F16_Nguyen_trim::trim_plant(vehicle.get_plant(), trim_inflags, trim_inputs); 

    auto [q,u] = trim2state(trim_output, trim_inputs);

    if ( !trim_success ) throw std::runtime_error("Trim did not succeed");

    // (5) Check trim point: transform to cppad + run operator()
    std::array<CppAD::AD<scalar>, F16_Nguyen::F16_Nguyen_plant::num_states()> q_cppad;
    std::array<CppAD::AD<scalar>, F16_optimal_control<>::NCONTROL> u_cppad;

    std::copy(q.cbegin(), q.cend(), q_cppad.begin());
    std::copy(u.cbegin(), u.cend(), u_cppad.begin());
    const auto dqdt = vehicle(q_cppad,{},u_cppad,0.0).first;
    for (const auto& dqdt_i : dqdt)
    {
        assert(std::abs(dqdt_i) < 1.0e-13 );
    }


    std::cout << "[info] Trim successful:" << std::endl;
    std::cout << "[info]     -> state: " << q << std::endl;
    std::cout << "[info]     -> controls: " << u << std::endl;
    std::cout << "[info]     -> time derivative: " << dqdt << std::endl;

  
    // (6) Run optimal control in a ovaltrack line
    auto control_variables = typename Optimal_laptime<F16_optimal_control<Track_by_arcs>>::template Control_variables<>{};

    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_full_mesh(std::vector<scalar>(n+1,u[F16_Nguyen::F16_Nguyen_plant::input_names::dh_dmd_deg]), 1.0e-4);

    vehicle.dlef_dmd_deg_default = u[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg];
    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::dlef_dmd_deg] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_dont_optimize();

    vehicle.dsb_dmd_deg_default = u[F16_Nguyen::F16_Nguyen_plant::input_names::dsb_dmd_deg];
    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::dsb_dmd_deg] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_dont_optimize();

    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::da_dmd_deg] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_full_mesh(std::vector<scalar>(n+1,u[F16_Nguyen::F16_Nguyen_plant::input_names::da_dmd_deg]), 1.0e-3);

    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::dr_dmd_deg] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_full_mesh(std::vector<scalar>(n+1,u[F16_Nguyen::F16_Nguyen_plant::input_names::dr_dmd_deg]), 1.0e-3);

    vehicle.thrustvectoring_longitude_angle_dmd_deg_default = u[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_longitude_angle_dmd_deg];
    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_longitude_angle_dmd_deg] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_dont_optimize();

    vehicle.thrustvectoring_latitude_angle_dmd_deg_default = u[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_latitude_angle_dmd_deg];
    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::thrustvectoring_latitude_angle_dmd_deg] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_dont_optimize();

    control_variables[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
        create_full_mesh(std::vector<scalar>(n+1,u[F16_Nguyen::F16_Nguyen_plant::input_names::throttle_percent]), 1.0e-4);

  
/*
    for (size_t i = 0; i < F16_optimal_control<>::NCONTROL; ++i)
        control_variables[i] = Optimal_laptime<F16_optimal_control<Track_by_arcs>>::
            create_full_mesh(std::vector<scalar>(n+1,u[i]), 1.0e-2);
*/
    

    std::vector<std::array<scalar,F16_optimal_control<>::NSTATE>> q_full(n+1,q);

    for (size_t i = 0; i < n+1; ++i)
    {
        const auto& pitch_deg = trim_output[F16_Nguyen::F16_Nguyen_trim::output_names::pitch_deg];
        const auto psi = vehicle.track_heading_angle_at(s[i]); 
        q_full[i][F16_Nguyen::F16_Nguyen_plant::state_names::qw_body2Earth] = cos(psi*0.5)*cos(pitch_deg*DEG*0.5);
        q_full[i][F16_Nguyen::F16_Nguyen_plant::state_names::qx_body2Earth] = -sin(psi*0.5)*sin(pitch_deg*DEG*0.5);
        q_full[i][F16_Nguyen::F16_Nguyen_plant::state_names::qy_body2Earth] = cos(psi*0.5)*sin(pitch_deg*DEG*0.5);
        q_full[i][F16_Nguyen::F16_Nguyen_plant::state_names::qz_body2Earth] = cos(pitch_deg*DEG*0.5)*sin(psi*0.5);
    }

    Optimal_laptime<F16_optimal_control<Track_by_arcs>>::Options opts;
    opts.print_level = 5;
    opts.maximum_iterations = 500;
    opts.throw_if_fail = false;
    opts.sigma = 0.5;
    Optimal_laptime opt_laptime(s, false, true, vehicle, q_full, {n+1,std::array<double,0>{}}, control_variables, opts);

    opt_laptime.xml()->save("results.xml");
}
