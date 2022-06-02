%
% Test an S & L trim, against the "boom3d" equivalent.
%

addpath(genpath('./functions'));


% open lib
libalias = 'libF16_Nguyen_clib';

dtor = onCleanup(@()(dlclose_F16_Nguyen_clib(libalias)));

dlreset_F16_Nguyen_clib(libalias, ...
    ['../../lib/Release/', libalias], ...
    '../../include/F16_Nguyen/F16_Nguyen_clib.h');


% create the plant
plant_properties = F16_Nguyen_plant_properties(libalias);

cplant = new_F16_Nguyen_plant(libalias, ...
    '../../datasets/aero_betasym', '../../datasets/engine');


% S & L trim definition
KCAS        = 250;
ZP_ft       = 15000;
gamma_deg   = 25;
disable_lef = false;
mass_t      = 9;
xcg_per_MAC = 0.55;


% "boom3d" solution
[boom3d_solution.theta_deg, ...
    boom3d_solution.dh_deg] = F16_pitch_for_trimmed_flight(KCAS, ZP_ft, ...
    0., gamma_deg, mass_t, xcg_per_MAC, disable_lef);

boom3d_solution.alpha_deg = boom3d_solution.theta_deg - gamma_deg;


% "F16_Nguyen_plant" solution
trim_inflags                         = plant_properties.default_trim_inflags;
trim_inflags.steady_trim             = true;
trim_inflags.straight_and_level_trim = true;
trim_inflags.disable_lef             = disable_lef;

trim_inputs                       = plant_properties.default_trim_inputs;
trim_inputs.KCAS                  = KCAS;
trim_inputs.ZP_ft                 = ZP_ft;
trim_inputs.flight_path_angle_deg = gamma_deg;
trim_inputs.mass_kg               = 1e3 * mass_t;
trim_inputs.xcg_per_MAC           = xcg_per_MAC;

[trim_success, ...
    trim_outputs] = F16_Nguyen_trim_plant(trim_inflags, trim_inputs, ...
    libalias, cplant, plant_properties);


% a solution using the Matlab function we have (which
% calls the plant's aero and engine datasets)
[matlab_solution.trim_success, ...
    matlab_solution.theta_deg, matlab_solution.dh_deg, ...
    matlab_solution.aoa_deg, matlab_solution.dlef_deg, ...
    matlab_solution.P3_percent, matlab_solution.throttle_percent, ...
    matlab_solution.costs] = F16_Nguyen_straight_and_level_trim(libalias, ...
    cplant, plant_properties, ...
    KCAS, ZP_ft, ...
    gamma_deg, mass_t, xcg_per_MAC, disable_lef);


% the two previous solutions should be super similar...
tol = 1e-10;
assert(isequal(trim_success, matlab_solution.trim_success))
assert(abs(trim_outputs.pitch_deg - matlab_solution.theta_deg) < tol)
assert(abs(trim_outputs.dh_deg - matlab_solution.dh_deg) < tol)
assert(abs(trim_outputs.aoa_deg - matlab_solution.aoa_deg) < tol)
assert(abs(trim_outputs.dlef_deg - matlab_solution.dlef_deg) < tol)
assert(abs(trim_outputs.P3_percent - matlab_solution.P3_percent) < tol)
assert(abs(trim_outputs.throttle_percent - matlab_solution.throttle_percent) < tol)


% summary
fprintf('test_straight_and_level_trim: summary:\n\n')
fprintf('  * F16_Nguyen trim success is %d.\n\n', trim_success)
fprintf('  * abs(F16_Nguyen_theta_deg - boom3d_theta_deg) is %.17g.\n\n', abs(trim_outputs.pitch_deg - boom3d_solution.theta_deg))
fprintf('  * abs(F16_Nguyen_dh_deg - boom3d_dh_deg) is %.17g.\n\n', abs(trim_outputs.dh_deg - boom3d_solution.dh_deg))
fprintf('  * F16_Nguyen''s costs are [%.17g, %.17g, %.17g].\n\n', [trim_outputs.cost_Fz, trim_outputs.cost_My, trim_outputs.cost_Fx])
fprintf('  * F16_Nguyen solution is:\n');
fprintf('    -> theta_deg        = %.17g.\n', trim_outputs.pitch_deg);
fprintf('    -> dh_deg           = %.17g.\n', trim_outputs.dh_deg);
fprintf('    -> aoa_deg          = %.17g.\n', trim_outputs.aoa_deg);
fprintf('    -> dlef_deg         = %.17g.\n', trim_outputs.dlef_deg);
fprintf('    -> P3_percent       = %.17g.\n', trim_outputs.P3_percent);
fprintf('    -> throttle_percent = %.17g.\n', trim_outputs.throttle_percent);


% finalize the plant and the lib
delete_F16_Nguyen_plant(libalias, cplant);
clear dtor