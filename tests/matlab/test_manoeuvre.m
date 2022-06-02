%
% Tests a manoeuvre.
%

addpath('./functions');
addpath(genpath('/home/c83833/projects/arbs/boom3d/visualize/functions'));


% helper constants
ft2m    = 0.3048;
kn2mps  = 0.514444;
deg2rad = pi / 180.;
rad2deg = 1. / deg2rad;


% open lib
libalias = 'libF16_Nguyen_clib';

dtor = onCleanup(@()(dlclose_F16_Nguyen_clib(libalias)));

dlreset_F16_Nguyen_clib(libalias, ...
    ['../../lib/Release/', libalias], ...
    '../../include/F16_Nguyen/F16_Nguyen_clib.h');


% create the plant
plant_properties = F16_Nguyen_plant_properties(libalias);

cplant = new_F16_Nguyen_plant(libalias, ...
    '../../datasets/aero', '../../datasets/engine');


% trim the plant
trim_inflags                         = plant_properties.default_trim_inflags;
trim_inflags.steady_trim             = true;
trim_inflags.straight_and_level_trim = true;
trim_inflags.disable_lef             = false;

trim_inputs       = plant_properties.default_trim_inputs;
trim_inputs.KCAS  = 270;
trim_inputs.ZP_ft = 27000;


[trim_success, ...
    trim_outputs] = F16_Nguyen_trim_plant(trim_inflags, trim_inputs, ...
    libalias, cplant, plant_properties);

if ~trim_success
    error('test_manoeuvre: unsuccessful trim...')
end


% set the initial states
x = cell2struct(num2cell(zeros(plant_properties.num_states, 1)), ...
    plant_properties.state_names);

ZP_m     = ft2m * trim_inputs.ZP_ft;
TAS_mps  = cas2tas(kn2mps * trim_inputs.KCAS, ZP_m);
x.u_mps = TAS_mps * cosd(trim_outputs.aoa_deg);
x.w_mps = TAS_mps * sind(trim_outputs.aoa_deg);

quat_body2Earth = ea2quat(deg2rad * [0; trim_outputs.pitch_deg; 0]);
x.qw_body2Earth = quat_body2Earth(1);
x.qx_body2Earth = quat_body2Earth(2);
x.qy_body2Earth = quat_body2Earth(3);
x.qz_body2Earth = quat_body2Earth(4);

x.zEarth_m   = -ZP_m;
x.dh_deg     = trim_outputs.dh_deg;
x.dlef_deg   = trim_outputs.dlef_deg;
x.P3_percent = trim_outputs.P3_percent;
x            = cell2mat(struct2cell(x));


% set the initial inputs
u                  = plant_properties.default_inputs;
u.dh_dmd_deg       = trim_outputs.dh_deg;
u.dlef_dmd_deg     = trim_outputs.dlef_deg;
u.throttle_percent = trim_outputs.throttle_percent;
u.mass_kg          = trim_inputs.mass_kg;
u.xcg_per_MAC      = trim_inputs.xcg_per_MAC;
u                  = cell2mat(struct2cell(u));


% calculate the initial outputs
y = F16_Nguyen_plant_outputs(x, u, ...
    libalias, cplant, plant_properties);


% simulate
sample_time_s = 0.016;
stop_time_s   = 1000.;
fcl           = struct(...
    'DNwz_dmd_g', 0., ...
    'Kp_long', -10.0, ...
    'Ki_long', -2, ...
    'integrator_long', trim_outputs.dh_deg, ...
    ...
    'beta_dmd_deg', 0., ...
    'Kp_dir', 0.5, ...
    'Ki_dir', 0.02, ...
    'integrator_dir', 0., ...
    ...
    'Dpw_dmd_degps', 0., ...
    'Kp_lat', -0.5, ...
    'Ki_lat', -0.02, ...
    'integrator_lat', 0., ...
    ...
    'TAS_dmd_mps', TAS_mps, ...
    'Kp_throttle', 5, ...
    'Ki_throttle', 0.5, ...
    'integrator_throttle', trim_outputs.throttle_percent);

num_ticks = stop_time_s / sample_time_s;

timehistories                       = struct();
timehistories.time_s                = zeros(num_ticks, 1);
timehistories.outputs(num_ticks, 1) = cell2struct(num2cell(zeros(plant_properties.num_outputs, 1)), ...
    plant_properties.output_names);
timehistories.inputs(num_ticks, 1)  = cell2struct(num2cell(zeros(plant_properties.num_inputs, 1)), ...
    plant_properties.input_names);
timehistories.states(num_ticks, 1)  = cell2struct(num2cell(zeros(plant_properties.num_states, 1)), ...
    plant_properties.state_names);

for tick = 1:num_ticks
    t = (tick - 1) * sample_time_s;

    [x, y] = continuous_plant_step(x, y, u, sample_time_s, ...
        libalias, cplant, plant_properties);
    
    [u, fcl] = ...
        update_plant_inputs(x, y, u, t, sample_time_s, fcl, ...
        plant_properties);
    
    timehistories.time_s(tick)     = t;
    timehistories.outputs(tick, 1) = cell2struct(num2cell(y), plant_properties.output_names);
    timehistories.inputs(tick, 1)  = cell2struct(num2cell(u), plant_properties.input_names);
    timehistories.states(tick, 1)  = cell2struct(num2cell(x), plant_properties.state_names);
    
end


% plot
figure
hold on
title('TAS mps')
plot(timehistories.time_s, [timehistories.outputs.TAS_mps])

figure
hold on
title('alpha deg')
plot(timehistories.time_s, [timehistories.outputs.aoa_deg])

figure
hold on
title('beta deg')
plot(timehistories.time_s, [timehistories.outputs.aos_deg])

figure
hold on
title('roll deg')
plot(timehistories.time_s, [timehistories.outputs.roll_deg])

figure
hold on
title('p degps')
plot(timehistories.time_s, rad2deg * [timehistories.states.p_radps])

figure
hold on
title('da deg')
plot(timehistories.time_s, [timehistories.states.da_deg])

figure
title('trajectory m')
hold on
axis equal
view(3)
xlabel('x')
ylabel('y')
zlabel('z')
plot3([timehistories.states.xEarth_m], ...
    [timehistories.states.yEarth_m], ...
    [timehistories.states.zEarth_m])


%
% Helper functions.
%

function [u, fcl] = update_plant_inputs(x, y, u, t, ...
    sample_time_s, fcl, ...
    plant_properties)
%
% Updates the plant inputs at each discrete step.
%

% helper constants
rad2deg = 180 / pi;


% outputs & inputs as struct
x = cell2struct(num2cell(x), plant_properties.state_names); %#ok<NASGU>
y = cell2struct(num2cell(y), plant_properties.output_names);
u = cell2struct(num2cell(u), plant_properties.input_names);


% longitudinal control
DNwz_g              = y.Nwz_g - cosd(y.flight_path_angle_deg) / cosd(y.flight_roll_angle_deg);
DNwz_error_g        = fcl.DNwz_dmd_g - DNwz_g;
fcl.integrator_long = fcl.integrator_long + fcl.Ki_long * sample_time_s * DNwz_error_g;
u.dh_dmd_deg        = fcl.Kp_long * DNwz_error_g + fcl.integrator_long;


% directional control
beta_error_deg     = fcl.beta_dmd_deg - y.aos_deg;
fcl.integrator_dir = fcl.integrator_dir + fcl.Ki_dir * sample_time_s * beta_error_deg;
u.dr_dmd_deg       = fcl.Kp_dir * beta_error_deg + fcl.integrator_dir;


% throttle control
TAS_error_mps           = fcl.TAS_dmd_mps - y.TAS_mps;
fcl.integrator_throttle = fcl.integrator_throttle + fcl.Ki_throttle * sample_time_s * TAS_error_mps;
u.throttle_percent      = fcl.Kp_throttle * TAS_error_mps + fcl.integrator_throttle;


% lateral input and control
if t > 1 && t <= 1.3
    u.da_dmd_deg = 21.5;

else
    Dpw_degps          = rad2deg * y.pw_radps + y.true_heading_angledot_degps * sind(y.flight_path_angle_deg);
    pw_error_degps     = fcl.Dpw_dmd_degps - Dpw_degps;
    fcl.integrator_lat = fcl.integrator_lat + fcl.Ki_lat * sample_time_s * pw_error_degps;
    u.da_dmd_deg       = fcl.Kp_lat * pw_error_degps + fcl.integrator_lat;
    
end

% inputs as a flat vector of doubles
u = cell2mat(struct2cell(u));

end


function [x, y] = continuous_plant_step(x, y, u, sample_time, ...
    libalias, cplant, plant_properties)
%
% Takes a continuous step of the plant
% with an RK4 scheme.
%

x0 = x;

k1 = F16_Nguyen_plant_derivatives(x, y, u, ...
    libalias, cplant, plant_properties);

[x, y] = update_plant_states(x0 + 0.5 * sample_time * k1, u, ...
    libalias, cplant, plant_properties);

k2 = F16_Nguyen_plant_derivatives(x, y, u, ...
    libalias, cplant, plant_properties);

[x, y] = update_plant_states(x0 + 0.5 * sample_time * k2, u, ...
    libalias, cplant, plant_properties);

k3 = F16_Nguyen_plant_derivatives(x, y, u, ...
    libalias, cplant, plant_properties);

[x, y] = update_plant_states(x0 + sample_time * k3, u, ...
    libalias, cplant, plant_properties);

k4 = F16_Nguyen_plant_derivatives(x, y, u, ...
    libalias, cplant, plant_properties);

[x, y] = ...
    update_plant_states(x0 + sample_time * (k1 + 2 * (k2 + k3) + k4) / 6, u, ...
    libalias, cplant, plant_properties);

end


function [x, y] = update_plant_states(x, u, ...
    libalias, cplant, plant_properties)
%
% Updtaes the states of the plant,
% returns the new states and the
% corresponding outputs.
%

y = F16_Nguyen_plant_outputs(x, u, ...
    libalias, cplant, plant_properties);

x = F16_Nguyen_plant_state_limiters(x, y, u, ...
    libalias, cplant, plant_properties);

y = F16_Nguyen_plant_outputs(x, u, ...
    libalias, cplant, plant_properties);

end