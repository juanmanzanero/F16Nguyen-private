function [trim_success, ...
    theta_deg, dh_deg, aoa_deg, dlef_deg, ...
    P3_percent, throttle_percent, ...
    trim_costs] = F16_Nguyen_straight_and_level_trim(libalias, cplant, plant_properties, ...
    KCAS, ZP_ft, gamma_deg, mass_t, xcg_per_MAC, disable_lef)
%
% Calculates the trimmed flight pitch of the "F16_Nguyen_plant", by solving
% the simplified longitudinal balance that assumes only two control
% surfaces (horizontal stabilator and leading edge flap):
%
%   q * S * CZ(aoa, dh, dlef) / W + cos(theta)      = 0,
%   CM(aoa, dh, dlef)                               = 0,
%   aoa + gamma - theta                             = 0,
%   dlef_deg - 1.38 * aoa_deg + 9.05 * q / p - 1.45 = 0,
%
% After balancing the aerodynamics, we can extract the engine's state
% "P3" and input "throttle" by solving:
%
%   Thrust(P3) / W - sin(theta) + q * S * CX(aoa, dh, dlef) / W = 0,
%   P3dot(throttle)                                             = 0.
%
% NOTE: this trim only evaluates the plant's aero-force and
% moment coefficients, not its state-space functions, and considers
% an engine totally aligned with the x-body-axis.
%
% Ref: "Simulator study of stall/post-stall characteristics of
% a fighter airplane with relaxed longitudinal static stability",
% L. T. Nguyen, NASA Technical Report, 1979.
%

% calculate the trim constants
[basic_aeroforce_and_moment_inputs, ...
    basic_engine_inputs, ...
    qS_div_W, q_div_p, W_N, ...
    theta_seed_deg, ...
    dh_seed_deg] = straight_and_level_flight_condition(libalias, cplant, plant_properties, ...
    KCAS, ZP_ft, gamma_deg, mass_t, xcg_per_MAC, disable_lef);


% set up the root finding algorithm for the aero-trim
cost_function = @(x)(straight_and_level_costs(x(1), x(2), ...
    libalias, cplant, plant_properties, ...
    gamma_deg, disable_lef, ...
    basic_aeroforce_and_moment_inputs, ...
    qS_div_W, q_div_p));

x0            = [theta_seed_deg; dh_seed_deg];
delta_numjac  = 0.1 * sqrt(eps);
alpha_numjac  = 1.;
max_fun_evals = 5000;
tolx          = 1e-16;
tolf          = 1e-14;

[x, costs_aero, exitflag] = numjac_Newton_method(cost_function, ...
    x0, ...
    delta_numjac, alpha_numjac, max_fun_evals, tolx, tolf);

aero_trim_success = exitflag == 1;
theta_deg         = x(1);
dh_deg            = x(2);

if ~aero_trim_success
    warning(['F16_Nguyen_straight_and_level_trim: the AERO root finding algorithm', ...
        ' failed, result will be approximate.'])

    seed_costs = cost_function(x0);
    if sumabs(seed_costs) < sumabs(costs_aero)
        costs_aero = seed_costs;
        theta_deg  = theta_seed_deg;
        dh_deg     = dh_seed_deg;

    end

end


% final aero outputs
aoa_deg  = theta_deg - gamma_deg;
dlef_deg = stationary_leading_edge_flap(aoa_deg, q_div_p, disable_lef);


% we now solve for the engine's state "P3" & input "throttle"
necessary_thrust_N = straight_and_level_thrust(theta_deg, ...
    aoa_deg, dh_deg, dlef_deg, ...
    libalias, cplant, plant_properties, ...
    basic_aeroforce_and_moment_inputs, ...
    qS_div_W, W_N);

[engine_trim_success, ...
    P3_percent, throttle_percent, ...
    cost_engine] = F16_Nguyen_engine_trim(libalias, cplant, plant_properties, ...
    basic_engine_inputs, ...
    necessary_thrust_N);


% finalize outputs
trim_success = aero_trim_success & engine_trim_success;
trim_costs   = [costs_aero(:); cost_engine(:)];

end


%
% Helper functions.
%

function costs = straight_and_level_costs(theta_deg, dh_deg, ...
    libalias, cplant, plant_properties, ...
    gamma_deg, disable_lef, ...
    basic_aeroforce_and_moment_inputs, ...
    qS_div_W, q_div_p)
%
% Longitudinal aerodynamic cost function for S & L flight.
%

basic_aeroforce_and_moment_inputs.aoa_deg = theta_deg - gamma_deg;
basic_aeroforce_and_moment_inputs.dh_deg  = dh_deg;

basic_aeroforce_and_moment_inputs.dlef_deg = ...
    stationary_leading_edge_flap(basic_aeroforce_and_moment_inputs.aoa_deg, ...
    q_div_p, disable_lef);

CFM = F16_Nguyen_aeroforce_and_moment_coefficients(basic_aeroforce_and_moment_inputs, ...
    libalias, cplant, plant_properties);

CZtot = CFM(strcmp('CZtot', plant_properties.aeroforce_and_moment_coefficient_names));
CMtot = CFM(strcmp('CMtot', plant_properties.aeroforce_and_moment_coefficient_names));

costs = [qS_div_W * CZtot + cosd(theta_deg); CMtot];

end


function dlef_deg = stationary_leading_edge_flap(aoa_deg, q_div_p, disable_lef)
%
% Stationary leading edge flap deflection, from the
% "Nguyen" reference.
%

if ~disable_lef
    dlef_deg = 1.38 * aoa_deg - 9.05 * q_div_p + 1.45;
else
    dlef_deg = 0.;
end

end


function [engine_trim_success, ...
    P3_percent, throttle_percent, ...
    cost_engine] = F16_Nguyen_engine_trim(libalias, cplant, plant_properties, ...
    basic_engine_inputs, ...
    necessary_thrust_N)
%
% x-body engine trim for S & L flight. The stationary solution
% "P3dot = 0" can only be achieved if "P1(throttle) = P2 = P3(thrust)",
% so we can obtain state "P3" and input "throttle" from the thrust.
%

% initialize output flag
engine_trim_success = true;


% calculate the engine's regime and "P3 = P3(thrust)"
engine_dataset_outputs = F16_Nguyen_engine_dataset(basic_engine_inputs, ...
    libalias, cplant, plant_properties);

Tidle_N = engine_dataset_outputs(strcmp('Tidle_N', plant_properties.engine_dataset_output_names));
Tmil_N  = engine_dataset_outputs(strcmp('Tmil_N', plant_properties.engine_dataset_output_names));
Tmax_N  = engine_dataset_outputs(strcmp('Tmax_N', plant_properties.engine_dataset_output_names));

P3_high_nondim = (necessary_thrust_N - Tmil_N) / (Tmax_N - Tmil_N) + 1.;
high_regime    = P3_high_nondim >= 1;

P3_low_nondim = (necessary_thrust_N - Tidle_N) / (Tmil_N - Tidle_N);
low_regime    = P3_low_nondim <= 1.;

if (high_regime && ~low_regime) || ...
        (high_regime && low_regime && P3_high_nondim == P3_low_nondim)
    
    P3_percent = P3_high_nondim * 50.;

elseif low_regime && ~high_regime
    P3_percent = P3_low_nondim * 50.;

else
    warning(['F16_Nguyen_straight_and_level_trim::F16_Nguyen_engine_trim:', ...
        ' the engine''s P3 regime is UNFEASIBLE, results will be approximate.']);

    engine_trim_success = false;
    
    in_high = P3_high_nondim - 1;
    in_low  = 1. - P3_low_nondim;
    if in_high > in_low
        P3_percent = P3_high_nondim * 50.;
        
    else
        P3_percent = P3_low_nondim * 50.;
        
    end

end

if P3_percent > 100. || P3_percent < 0.
    warning(['F16_Nguyen_straight_and_level_trim::F16_Nguyen_engine_trim:', ...
        ' the engine''s P3 regime is UNFEASIBLE, results will be approximate.']);

    engine_trim_success = false;
    P3_percent          = max(0., min(P3_percent, 100));
    throttle_percent    = P3_percent;
    cost_engine         = inf;
    return

end


% extract the throttle from "P3 = P1(throttle)",
% solving the curve with fzero in the interval
% "[0, 100]"
tolf          = 1e-14;
cost_function = @(x)(engine_throttle_cost(x, ...
    libalias, cplant, plant_properties, ...
    basic_engine_inputs, ...
    P3_percent));

[x, cost_engine, exitflag] = ...
    fzero(cost_function, [0 100], struct('TolX', tolf));

if exitflag ~= 1
    warning(['F16_Nguyen_straight_and_level_trim::F16_Nguyen_engine_trim:', ...
        ' the THROTTLE root finding algorithm failed, result will be approximate.'])

    engine_trim_success = false;
    throttle_percent    = P3_percent;

else
    throttle_percent = x;

end

end


function cost = engine_throttle_cost(throttle_percent, ...
    libalias, cplant, plant_properties, ...
    basic_engine_inputs, ...
    necessary_P1_percent)
%
% Cost function to solve the engine's input throttle
% from a known "P1_percent".
%

basic_engine_inputs.throttle_percent = throttle_percent;

engine_dataset_outputs = F16_Nguyen_engine_dataset(basic_engine_inputs, ...
    libalias, cplant, plant_properties);

P1_percent = engine_dataset_outputs(strcmp('P1_percent', ...
    plant_properties.engine_dataset_output_names));

cost = P1_percent - necessary_P1_percent;

end


function [basic_aeroforce_and_moment_inputs, ...
    basic_engine_inputs, ...
    qS_div_W, q_div_p, W_N, ...
    theta_point_deg, ...
    dh_point_deg] = straight_and_level_flight_condition(libalias, cplant, plant_properties, ...
    KCAS, ZP_ft, gamma_deg, mass_t, xcg_per_MAC, disable_lef)
%
% Calculates some flight-condition-related constants
% for the S & L flight.
%

% helper constants
ft2m    = 0.3048;
kn2mps  = 0.514444;
t2kg    = 1e3;
g0_mps2 = 9.80665;


% basic inputs to the plant aero model,
% basic inputs & states for the engine model
basic_aeroforce_and_moment_inputs = ...
    cell2struct(num2cell(zeros(plant_properties.num_aeroforce_and_moment_inputs, 1)), ...
    plant_properties.aeroforce_and_moment_input_names);

basic_engine_inputs = ...
    cell2struct(num2cell(zeros(plant_properties.num_engine_inputs, 1)), ...
    plant_properties.engine_input_names);

for di = fieldnames(plant_properties.default_inputs)'
    if ismember(di{:}, fieldnames(basic_aeroforce_and_moment_inputs))
        basic_aeroforce_and_moment_inputs.(di{:}) = ...
            plant_properties.default_inputs.(di{:});

    end

    if ismember(di{:}, fieldnames(basic_engine_inputs))
        basic_engine_inputs.(di{:}) = ...
            plant_properties.default_inputs.(di{:});

    end

end


% calculate flight condition
basic_engine_inputs.ZP_m                           = ft2m * ZP_ft;
[density_kgpm3, soundspeed_mps, ~, ~, pressure_Pa] = isa_atmos(basic_engine_inputs.ZP_m, 0.);
basic_aeroforce_and_moment_inputs.xcg_per_MAC      = xcg_per_MAC;
basic_aeroforce_and_moment_inputs.TAS_mps          = cas2tas(kn2mps * KCAS, basic_engine_inputs.ZP_m, 0.);
basic_engine_inputs.Mach                           = basic_aeroforce_and_moment_inputs.TAS_mps / soundspeed_mps;
qbar_Pa                                            = 0.5 * density_kgpm3 * basic_aeroforce_and_moment_inputs.TAS_mps^2;
W_N                                                = g0_mps2 * t2kg * mass_t;
qS_div_W                                           = qbar_Pa * plant_properties.default_inputs.wing_surface_m2 / W_N;
q_div_p                                            = qbar_Pa / pressure_Pa;


% calculate the theta seed with the simplest
% wind-axes point formula
seed_aeroforce_and_moment_inputs = basic_aeroforce_and_moment_inputs;

Daoa_deg      = 0.1 * sqrt(eps);
aoa_point_deg = 0.;
for iter_point = 1:2
    seed_aeroforce_and_moment_inputs.aoa_deg  = aoa_point_deg;

    seed_aeroforce_and_moment_inputs.dlef_deg = ...
        stationary_leading_edge_flap(seed_aeroforce_and_moment_inputs.aoa_deg, ...
        q_div_p, disable_lef);

    approx_CL = F16_Nguyen_total_lift_coefficient(seed_aeroforce_and_moment_inputs, ...
            libalias, cplant, plant_properties);

    approx_CLalpha = 0.;
    for s = [-1, 1]
        seed_aeroforce_and_moment_inputs.aoa_deg  = aoa_point_deg + s * Daoa_deg;

        seed_aeroforce_and_moment_inputs.dlef_deg = ...
            stationary_leading_edge_flap(seed_aeroforce_and_moment_inputs.aoa_deg, ...
            q_div_p, disable_lef);

        approx_CLalpha = approx_CLalpha + ...
            s *  F16_Nguyen_total_lift_coefficient(seed_aeroforce_and_moment_inputs, ...
            libalias, cplant, plant_properties);
        
    end

    approx_CLalpha = 0.5 * approx_CLalpha / Daoa_deg;
    approx_CL0     = approx_CL - approx_CLalpha * aoa_point_deg;

    aoa_point_deg = (cosd(gamma_deg) / qS_div_W - approx_CL0) / approx_CLalpha;

end

theta_point_deg = aoa_point_deg + gamma_deg;


% calculate the stabilator seed with the simplest
% longitudinal formula, 0 = CM0(alpha, lef, dh0) + CMdh * dh
seed_aeroforce_and_moment_inputs = basic_aeroforce_and_moment_inputs;

seed_aeroforce_and_moment_inputs.aoa_deg = aoa_point_deg;

seed_aeroforce_and_moment_inputs.dlef_deg = ...
    stationary_leading_edge_flap(seed_aeroforce_and_moment_inputs.aoa_deg, ...
    q_div_p, disable_lef);
    
Dh_deg        = 0.1 * sqrt(eps);
dh_point_deg  = 0.;
for iter_point = 1:2
    seed_aeroforce_and_moment_inputs.dh_deg = dh_point_deg;

    approx_CM = F16_Nguyen_total_pitch_moment_coefficient(seed_aeroforce_and_moment_inputs, ...
        libalias, cplant, plant_properties);
    
    approx_CMdh = 0.;
    for s = [-1, 1]
        seed_aeroforce_and_moment_inputs.dh_deg = dh_point_deg + s * Dh_deg;
        
        approx_CMdh = approx_CMdh + ...
            s * F16_Nguyen_total_pitch_moment_coefficient(seed_aeroforce_and_moment_inputs, ...
            libalias, cplant, plant_properties);
        
    end

    approx_CMdh = 0.5 * approx_CMdh / Dh_deg;
    approx_CM0  = approx_CM - approx_CMdh * dh_point_deg;

    dh_point_deg = -approx_CM0 / approx_CMdh;

end

end


function necessary_thrust_N = straight_and_level_thrust(theta_deg, ...
    aoa_deg, dh_deg, dlef_deg, ...
    libalias, cplant, plant_properties, ...
    basic_aeroforce_and_moment_inputs, ...
    qS_div_W, W_N)
%
% Calculates the necessary "T / W" for
% a S & L flight.
%

basic_aeroforce_and_moment_inputs.aoa_deg  = aoa_deg;
basic_aeroforce_and_moment_inputs.dh_deg   = dh_deg;
basic_aeroforce_and_moment_inputs.dlef_deg = dlef_deg;

CFM = F16_Nguyen_aeroforce_and_moment_coefficients(basic_aeroforce_and_moment_inputs, ...
    libalias, cplant, plant_properties);

CXtot = CFM(strcmp('CXtot', plant_properties.aeroforce_and_moment_coefficient_names));

necessary_thrust_N = W_N * (sind(theta_deg) - qS_div_W * CXtot);

end


function CFM = F16_Nguyen_aeroforce_and_moment_coefficients(aeroforce_and_moment_inputs, ...
    libalias, cplant, plant_properties)
%
% Calculates the force & moment coefficients for the given
% "aeroforce_and_moment_inputs".
%

[err, CFM] = calllib(libalias, 'F16_Nguyen_clib_plant_aeroforce_and_moment_coefficients', ...
    zeros(plant_properties.num_aeroforce_and_moment_coefficients, 1), plant_properties.num_aeroforce_and_moment_coefficients, ...
    cplant, ...
    cell2mat(struct2cell(aeroforce_and_moment_inputs)), ...
    plant_properties.num_aeroforce_and_moment_inputs);

if err ~= 0
    error(['F16_Nguyen_straight_and_level_trim::F16_Nguyen_aeroforce_and_moment_coefficients:', ...
        ' "F16_Nguyen_clib_plant_aeroforce_and_moment_coefficients" returned an', ...
        ' error (%d).'], err);

end

end


function CLtot = F16_Nguyen_total_lift_coefficient(aeroforce_and_moment_inputs, ...
    libalias, cplant, plant_properties)
%
% Calculates the lift coefficient of the plant for the given
% "aeroforce_and_moment_inputs".
%

CFM = F16_Nguyen_aeroforce_and_moment_coefficients(aeroforce_and_moment_inputs, ...
    libalias, cplant, plant_properties);

CZtot = CFM(strcmp('CZtot', plant_properties.aeroforce_and_moment_coefficient_names));
CXtot = CFM(strcmp('CXtot', plant_properties.aeroforce_and_moment_coefficient_names));
CLtot = -cosd(aeroforce_and_moment_inputs.aoa_deg) * CZtot + ...
    sind(aeroforce_and_moment_inputs.aoa_deg) * CXtot;

end


function CMtot = F16_Nguyen_total_pitch_moment_coefficient(aeroforce_and_moment_inputs, ...
    libalias, cplant, plant_properties)
%
% Calculates the pitch moment coefficient of the plant for the given
% "aeroforce_and_moment_inputs".
%

CFM = F16_Nguyen_aeroforce_and_moment_coefficients(aeroforce_and_moment_inputs, ...
    libalias, cplant, plant_properties);

CMtot = CFM(strcmp('CMtot', plant_properties.aeroforce_and_moment_coefficient_names));

end


function engine_dataset_outputs = F16_Nguyen_engine_dataset(engine_inputs, ...
    libalias, cplant, plant_properties)
%
% Evaluates the engine dataset of the "F16_Nguyen_plant",
% for the given "engine_inputs".
%

[err, engine_dataset_outputs] = calllib(libalias, 'F16_Nguyen_clib_plant_engine_dataset', ...
    zeros(plant_properties.num_engine_dataset_outputs, 1), plant_properties.num_engine_dataset_outputs, ...
    cplant, ...
    cell2mat(struct2cell(engine_inputs)), plant_properties.num_engine_inputs);

if err ~= 0
    error(['F16_Nguyen_straight_and_level_trim::F16_Nguyen_engine_dataset:', ...
        ' "F16_Nguyen_clib_plant_engine_dataset" returned an', ...
        ' error (%d).'], err);

end

end