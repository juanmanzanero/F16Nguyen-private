function steady_outputs = flat_Earth_steady_flight(steady_inputs, ...
    num_ticks, sample_time_s)
%
% Returns acceleration, velocity & attitude timehistories (and
% other useful constants that define the manoeuvre) for a
% steady flight (assuming the flat Earth model) with CONSTANT
% and SPECIFIED:
%
%   steady_inputs = struct('TAS_mps', , ... % true airspeed (scalar)
%       'alpha_deg', , ... % angle of attack
%       'beta_deg', , ... % angle of sideslip
%       'roll_deg', , ... % body roll angle
%       'gamma_deg' | 'pitch_deg', , ... % either the (airspeed) flight path angle or the body pitch angle
%       ...
%       'Ny_g' | 'Nz_g' | 'Nx_g' | 'Nwy_g' | 'Nwz_g', , ... % a load factor, it can be one of the three body-axes ones, or the "y" or
%       ...           % "z" wind-axes factor (optional, by default specifies "Nwy_g = 0", which is
%       ...           % the typical condition in steady turns). Specifying "Nwx_g" yields a
%       ...           % singular problem
%       ...
%       'initial_yaw_deg', , ... % the yaw angle at the first instant (optional, zero by default)
%       'wind_velocity_Earthaxes_mps', ); % mean wind velocity vector, projected onto the Earth-axes frame (optional, [0.; 0.; 0.] by default)
%
% This specification of inputs covers all of the classic
% steady flights:
%
%   * straight and level flight,
%   * steady turn,
%   * straight climb/descent,
%   * helicoidal flight (i.e., the turning climb/descent).
%

% helper constants
rad2deg = 180. / pi;
g0_mps2 = 9.80665;


% validate inputs
steady_inputs = validate_steady_inputs(steady_inputs);


% start by calculating some helper sines and cosines
calpha = cosd(steady_inputs.alpha_deg);
salpha = sind(steady_inputs.alpha_deg);
cbeta  = cosd(steady_inputs.beta_deg);
sbeta  = sind(steady_inputs.beta_deg);
cphi   = cosd(steady_inputs.roll_deg);
sphi   = sind(steady_inputs.roll_deg);

if isfield(steady_inputs, 'gamma_deg')
    % solve for the pitch, the problem has two
    % possible solutions (we'll take the one closest
    % to alpha), and it may be unfeasible for the
    % specified inputs (in which case we'll change
    % the roll angle)
    gamma_deg = steady_inputs.gamma_deg;
    cgamma    = cosd(steady_inputs.gamma_deg);
    
    [pitch_deg, feasible, ...
        roll_deg, cphi, sphi] = steady_flight_pitch(cgamma, sind(gamma_deg), ...
        steady_inputs.roll_deg, cphi, sphi, ...
        steady_inputs.alpha_deg, calpha, salpha, ...
        cbeta, sbeta);
    
    ctheta = cosd(pitch_deg);
    stheta = sind(pitch_deg);
    
else
    % solve for gamma, this problem is always well-posed
    feasible = true;
    
    roll_deg = steady_inputs.roll_deg;
    
    pitch_deg = steady_inputs.pitch_deg;
    ctheta    = cosd(pitch_deg);
    stheta    = sind(pitch_deg);
    
    gamma_deg = steady_flight_path_angle(ctheta, stheta, ...
        calpha, salpha, ...
        cbeta, sbeta, ...
        cphi, sphi);
    
    cgamma = cosd(gamma_deg);
    
end


% the wind-axes roll angle can be obtained by rotating
% the z-Earth vector to the wind-axes frame and looking
% at its second and third components
ctheta_cphi        = ctheta * cphi;
ctheta_sphi        = ctheta * sphi;
calpha_ctheta_cphi = calpha * ctheta_cphi;
calpha_stheta      = calpha * stheta;
salpha_ctheta_cphi = salpha * ctheta_cphi;
salpha_stheta      = salpha * stheta;
cgamma_smu         = cbeta * ctheta_sphi + sbeta * (calpha_stheta - salpha_ctheta_cphi);
cgamma_cmu         = calpha_ctheta_cphi + salpha_stheta;
mu_deg             = atan2d(cgamma_smu, cgamma_cmu);


% from alpha and beta we can obtain the true airspeed
% vector, projected in body-axes
VTAS_mps = steady_inputs.TAS_mps * [calpha * cbeta; ...
    sbeta;  salpha * cbeta];


% the yaw rate is also constant, and we can obtain it
% from the load factor equations. It can be shown that
% the wind velocity terms cancel out, because with
% (constant) wind, the body pseudo-acceleration
% isn't zero:
%
%   uvwvdot = -cross(pqr, wind_velocity_bodyaxes_mps).
%
% This condition means that in the body-load factor equations
% (flat Earth, g = [0; 0; g0]) we can take the "VTAS"
% directly (which is constant) instead of the inertial
% components "uvw", because their wind-related terms cancel
% out with "uvwdot". The wind-load factor equations are just
% the body-axes ones but rotated by "T_body2wind"
if isfield(steady_inputs, 'Nx_g')
    num_mps2 = g0_mps2 * (steady_inputs.Nx_g + stheta);
    den_mps  = VTAS_mps(2) * ctheta_cphi - VTAS_mps(3) * ctheta_sphi;

    yawdot_radps = ...
        check_yawdot_solution_with_specified_load_factor(num_mps2, den_mps, ...
        'Nx_g', '-sin(theta)');

elseif isfield(steady_inputs, 'Ny_g')
    num_mps2 = g0_mps2 * (ctheta_sphi - steady_inputs.Ny_g);
    den_mps  = VTAS_mps(1) * ctheta_cphi + VTAS_mps(3) * stheta;

    yawdot_radps = ...
        check_yawdot_solution_with_specified_load_factor(num_mps2, den_mps, ...
        'Ny_g', 'cos(theta) * sin(phi)');

elseif isfield(steady_inputs, 'Nz_g')
    num_mps2 = g0_mps2 * (steady_inputs.Nz_g - ctheta_cphi);
    den_mps  = VTAS_mps(1) * ctheta_sphi + VTAS_mps(2) * stheta;

    yawdot_radps = ...
        check_yawdot_solution_with_specified_load_factor(num_mps2, den_mps, ...
        'Nz_g', 'cos(theta) * cos(phi)');

elseif isfield(steady_inputs, 'Nwx_g')
    % this one is always singular (unless "Nwx_g == -sind(gamma_deg)")
    num_mps2 = g0_mps2 * (steady_inputs.Nwx_g + sind(gamma_deg));

    % den_zero_by_definition_mps = (cbeta * salpha * ctheta_sphi - sbeta * ctheta_cphi) * ...
    %     VTAS_mps(1) + cbeta * cgamma_cmu * VTAS_mps(2) - ...
    %     (sbeta * stheta + calpha * cbeta * ctheta_sphi) * VTAS_mps(3);
    den_zero_by_definition_mps = 0.;

    yawdot_radps = ...
        check_yawdot_solution_with_specified_load_factor(num_mps2, ...
        den_zero_by_definition_mps, ...
        'Nxw_g', '-sin(gamma)');

elseif isfield(steady_inputs, 'Nwy_g')
    num_mps2 = g0_mps2 * (cgamma_smu - steady_inputs.Nwy_g);
    den_mps  = steady_inputs.TAS_mps * (salpha_stheta + calpha_ctheta_cphi);

    yawdot_radps = ...
        check_yawdot_solution_with_specified_load_factor(num_mps2, den_mps, ...
        'Nwy_g', 'cos(gamma) * sin(mu)');

elseif isfield(steady_inputs, 'Nwz_g')
    num_mps2 = g0_mps2 * (steady_inputs.Nwz_g - cgamma_cmu);
    den_mps  = steady_inputs.TAS_mps * (sbeta * (calpha_stheta - ...
        salpha_ctheta_cphi) + cbeta * ctheta_sphi);

    yawdot_radps = ...
        check_yawdot_solution_with_specified_load_factor(num_mps2, den_mps, ...
        'Nwz_g', 'cos(gamma) * cos(mu)');

end


% angular velocity in body axes
z_Earth_bodyaxes = [-stheta; ctheta_sphi; ctheta_cphi];
pqr_radps        = yawdot_radps * z_Earth_bodyaxes;


% the inertial acceleration in body axes is constant, as said
% before, with wind, the body pseudo-accelerations "uvwdot"
% are nonzero and cancel out with wind terms of the Coriolis term
axyz_mps2 = cross(pqr_radps, VTAS_mps);


% the body and wind load factors are also constant, the latter
% are just the former but rotated by "T_body2wind"
Nxyz_g = z_Earth_bodyaxes - axyz_mps2 / g0_mps2;

T_body2wind = [calpha * cbeta, sbeta,  salpha * cbeta; ...
    -calpha * sbeta, cbeta, -salpha * sbeta; ...
    -salpha,     0,          calpha];

Nwxyz_g = T_body2wind * Nxyz_g;


% the "air" turn radius (i.e., the turn radius formed
% with the VTAS vector instead of the inertial speed)
% is constant, since chidot = yawdot (the angular velocity
% is "chidot * z_Earthaxes = yawdot * z_Earthaxes" because
% alpha and beta are constant).
air_turn_radius_m = steady_inputs.TAS_mps / ...
    abs(yawdot_radps * cgamma);


% timehistories: yaw is simply linear, and we'll generate
% speed & acceleration in Earth axes from it
yawdot_degps = rad2deg * yawdot_radps;

yaw_deg_timehistory = steady_inputs.initial_yaw_deg + ...
    yawdot_degps * sample_time_s * (0:num_ticks - 1)';

axyz_Earthaxes_mps2_timehistory   = zeros(num_ticks, 3);
uvwdot_mps2_timehistory           = zeros(num_ticks, 3);
uvw_Earthaxes_mps_timehistory     = zeros(num_ticks, 3);
uvw_mps_timehistory               = zeros(num_ticks, 3);
VTAS_Earthaxes_mps_timehistory    = zeros(num_ticks, 3);

for tick = 1:num_ticks
    [axyz_Earthaxes_mps2_timehistory(tick, :), ...
        uvwdot_mps2_timehistory(tick, :), ...
        uvw_Earthaxes_mps_timehistory(tick, :), ...
        uvw_mps_timehistory(tick, :), ...
        VTAS_Earthaxes_mps_timehistory(tick, :)] = ...
        steady_flight_velocity_and_acceleration(yaw_deg_timehistory(tick, 1), ...
        axyz_mps2, VTAS_mps, pqr_radps, pitch_deg, roll_deg, ...
        steady_inputs.wind_velocity_Earthaxes_mps);
    
end


% gather outputs
steady_inputs = rmfield(steady_inputs, 'roll_deg');
steady_inputs = try_rmfield(steady_inputs, 'gamma_deg');
steady_inputs = try_rmfield(steady_inputs, 'pitch_deg');
steady_inputs = try_rmfield(steady_inputs, 'Nx_g');
steady_inputs = try_rmfield(steady_inputs, 'Ny_g');
steady_inputs = try_rmfield(steady_inputs, 'Nz_g');

steady_outputs = merge_structs(steady_inputs, ...
    struct('steady_inputs_are_feasible', feasible, ...
    'gamma_deg', gamma_deg, ...
    'pitch_deg', pitch_deg, ...
    'roll_deg', roll_deg, ...
    'mu_deg', mu_deg, ...
    'VTAS_mps', VTAS_mps, ...
    'yawdot_degps', yawdot_degps, ...
    'pqr_degps', rad2deg * pqr_radps, ...
    'axyz_mps2', axyz_mps2, ...
    'Nxyz_g', Nxyz_g, ...
    'Nwxyz_g', Nwxyz_g, ...
    'air_turn_radius_m', air_turn_radius_m, ...
    ...
    'yaw_deg_timehistory', yaw_deg_timehistory, ...
    'axyz_Earthaxes_mps2_timehistory', axyz_Earthaxes_mps2_timehistory, ...
    'uvwdot_mps2_timehistory', uvwdot_mps2_timehistory, ...
    'uvw_Earthaxes_mps_timehistory', uvw_Earthaxes_mps_timehistory, ...
    'uvw_mps_timehistory', uvw_mps_timehistory, ...
    'VTAS_Earthaxes_mps_timehistory', VTAS_Earthaxes_mps_timehistory));

end


%
% Helper functions.
%

function steady_inputs = validate_steady_inputs(steady_inputs)
%
% Validates the structure of inputs to "flat_Earth_steady_flight".
%

% TAS
if steady_inputs.TAS_mps <= 0
    error(['flat_Earth_steady_flight::validate_inputs: "steady_inputs.TAS_mps"', ...
        ' should be positive.']);
    
end


% alpha
if abs(steady_inputs.alpha_deg) > 90
    error(['flat_Earth_steady_flight::validate_inputs: "steady_inputs.alpha_deg"', ...
        ' should be in [-90, 90].']);
    
end


% beta
steady_inputs.beta_deg = wrap_to_180(steady_inputs.beta_deg);


% body roll
steady_inputs.roll_deg = wrap_to_180(steady_inputs.roll_deg);


% gamma or pitch
if isfield(steady_inputs, 'gamma_deg') && ~isfield(steady_inputs, 'pitch_deg')
    if abs(steady_inputs.gamma_deg) > 90
        error(['flat_Earth_steady_flight::validate_inputs: "steady_inputs.gamma_deg"', ...
            ' should be in [-90, 90].']);
        
    end
    
elseif isfield(steady_inputs, 'pitch_deg') && ~isfield(steady_inputs, 'gamma_deg')
    if abs(steady_inputs.pitch_deg) > 90
        error(['flat_Earth_steady_flight::validate_inputs: "steady_inputs.pitch_deg"', ...
            ' should be in [-90, 90].']);
        
    end
    
else
    error(['flat_Earth_steady_flight::validate_inputs: either "steady_inputs.gamma_deg" or', ...
        ' "steady_inputs.pitch_deg" should be specified, and never both at once.']);
    
end


% load factors
specified_load_factors = isfield(steady_inputs, 'Nx_g') + ...
    isfield(steady_inputs, 'Ny_g') + isfield(steady_inputs, 'Nz_g') + ...
    isfield(steady_inputs, 'Nwx_g') + ...
    isfield(steady_inputs, 'Nwy_g') + isfield(steady_inputs, 'Nwz_g');

if specified_load_factors == 0
    steady_inputs.Nwy_g = 0.;
    
elseif specified_load_factors ~= 1
    error(['flat_Earth_steady_flight::validate_inputs: "steady_inputs" should', ...
        ' only specify one of the load factors (body or wind), i.e., either "Nx_g",', ...
        ' or "Ny_g", or "Nz_g", or "Nwx_g", or "Nwy_g", or "Nwz_g".']);
    
end


% initial yaw
if ~isfield(steady_inputs, 'initial_yaw_deg')
    steady_inputs.initial_yaw_deg = 0.;
else
    steady_inputs.initial_yaw_deg = wrap_to_180(steady_inputs.initial_yaw_deg);
end


% wind velocity
if ~isfield(steady_inputs, 'wind_velocity_Earthaxes_mps')
    steady_inputs.wind_velocity_Earthaxes_mps = [0; 0; 0];
    
else
    assert(numel(steady_inputs.wind_velocity_Earthaxes_mps) == 3)
    steady_inputs.wind_velocity_Earthaxes_mps = steady_inputs.wind_velocity_Earthaxes_mps(:);
    
end

end


function [pitch_deg, feasible, ...
    roll_deg, cphi, sphi] = steady_flight_pitch(cgamma, sgamma, ...
    roll_deg, cphi, sphi, ...
    alpha_deg, calpha, salpha, ...
    cbeta, sbeta)
%
% Returns the pitch of the steady flight when alpha, beta, roll
% and gamma are known. If the inputs are unfeasible, it
% re-calculates the roll (see below).
%

% helper constants
rad2deg = 180 / pi;


% solve the wind-angles relationship for the pitch
[thetas_rad, feasible] = sin_cos_solve(calpha * cbeta, ...
    -sphi * sbeta - cphi * salpha * cbeta, ...
    sgamma);

if feasible
    % since there are two PERFECTLY POSSIBLE solutions we'll follow
    % a criteria: have the pitch be similar to alpha, this seems the
    % most natural one...
    thetas_deg        = rad2deg * thetas_rad;
    [~, most_similar] = min(abs(thetas_deg - alpha_deg));
    pitch_deg         = thetas_deg(most_similar);
    
else
    % there are certain constraints that the specified [alpha,
    % beta, gamma, phi] must satisfy (see function "sin_cos_solve"),
    % so in fucked up situations we may end up here... for example
    % in turns of zero radius. In this case we'll do "mu <-- roll",
    % and recalculate the roll, to get some kind of solution...
    warning(['flat_Earth_steady_flight::steady_flight_pitch: inputs', ...
        ' are ill-posed, the specified "steady_inputs.roll_deg" will', ...
        ' be interpreted as "mu_deg" (wind axes roll) and the body-roll', ...
        ' angle will be re-calculated.']);
    
    cmu       = cphi;
    smu       = sphi;
    pitch_deg = asind(sgamma * calpha * cbeta + ...
        cgamma * smu * calpha * sbeta + cgamma * cmu * salpha);
    
    % unless we have gymbal lock, the roll is defined; else,
    % we'll arbitrarily take roll = 0.
    gymbal_lock_pitch_tol_deg = 1e-3;
    if abs(abs(pitch_deg) - 90.) > gymbal_lock_pitch_tol_deg
        roll_deg = wrap_to_180(atan2d(cgamma * smu * cbeta - sgamma * sbeta, ...
            cgamma * cmu * calpha - sgamma * salpha * cbeta - ...
            cgamma * smu * salpha * sbeta));
        
    else
        roll_deg = 0.;
        
    end
    
    % update the roll's cosine and sine
    cphi = cosd(roll_deg);
    sphi = sind(roll_deg);
    
end

end


function gamma_deg = steady_flight_path_angle(ctheta, stheta, ...
    calpha, salpha, ...
    cbeta, sbeta, ...
    cphi, sphi)
%
% Returns the (airspeed) flight path angle (gamma) of the
% steady flight when alpha, beta, roll and pitch are known.
%

% solve the wind-angles relationship for gamma
sgamma = calpha * cbeta * stheta - ...
    (sphi * sbeta + cphi * salpha * cbeta) * ctheta;

% clip it just for numerical reasons, the analytical result
% of the expression above is always in [-1, 1]
gamma_deg = asind(max(-1., min(1., sgamma)));

end


function [axyz_Earthaxes_mps2, ...
    uvwdot_mps2, ...
    uvw_Earthaxes_mps, ...
    uvw_mps, ...
    VTAS_Earthaxes_mps] = steady_flight_velocity_and_acceleration(yaw_deg, ...
    axyz_mps2, VTAS_mps, pqr_radps, pitch_deg, roll_deg, ...
    wind_velocity_Earthaxes_mps)
%
% Returns the inertial acceleration, pseudo-acceleration
% and velocity of the steady flight.
%

% helper constants
deg2rad = pi / 180.;


% inertial acceleration in Earth axes (we'll return a row vector,
% for convenience). Only changes due to the varying yaw.
T_body2Earth        = ea2rotmat(deg2rad * [yaw_deg; pitch_deg; roll_deg]);
axyz_Earthaxes_mps2 = (T_body2Earth * axyz_mps2)';


% pseudo-acceleration in body axes (we'll return a row vector,
% for convenience). Would be zero if there was no wind
wind_velocity_bodyaxes_mps = T_body2Earth' * wind_velocity_Earthaxes_mps;
uvwdot_mps2                = (-cross(pqr_radps, wind_velocity_bodyaxes_mps))';


% inertial velocity in body axes (we'll return a row vector,
% for convenience). Would be constant and equal to "VTAS_mps"
% if there was no wind
uvw_mps = (VTAS_mps + wind_velocity_bodyaxes_mps)';


% inertial velocity in Earth axes (we'll return a row vector,
% for convenience). Only changes due to the varying yaw
uvw_Earthaxes_mps = (T_body2Earth * uvw_mps')';


% true airspeed vector in Earth axes (we'll return a row vector,
% for convenience). Only changes due to the varying yaw
VTAS_Earthaxes_mps = (T_body2Earth * VTAS_mps)';

end


function yawdot_radps = ...
    check_yawdot_solution_with_specified_load_factor(num_mps2, den_mps, ...
    Nxyz_str, zero_num_kinematic_solution_str)
%
% Checks if a "yawdot_radps = num_mps2 / den_mps" solution with specified
% "N(w)z", "N(w)y",  or "N(w)z" is feasible, and returns its value.
%

tol = 1e3 * eps;

if abs(num_mps2) < tol
    if abs(den_mps) < tol
        warning(['flat_Earth_steady_flight::check_yawdot_solution_with_specified_load_factor:', ...
            ' specifying "', Nxyz_str,'" yields an ill-posed problem, but since', ...
            ' the input "', Nxyz_str,'" is equal to "', zero_num_kinematic_solution_str, '",', ...
            ' a solution with zero yaw rate exists.']);
        
    end
    
    yawdot_radps = 0.;
    
elseif abs(den_mps) < tol
    error(['flat_Earth_steady_flight::check_yawdot_solution_with_specified_load_factor:', ...
        ' specifying "', Nxyz_str, '" yields a singular problem.']);
    
else
    yawdot_radps = num_mps2 / den_mps;
    
end

end


function T = ea2rotmat(ea_rad)
%
% Rotation matrix from 'ZYX' (yaw-pitch-roll) Euler angle sequence [rad].
% The rotation matrix transforms from the rotated axes to the original
% ones: x_original = T * x_rotated.
%

cpsi   = cos(ea_rad(1));
spsi   = sin(ea_rad(1));
ctheta = cos(ea_rad(2));
stheta = sin(ea_rad(2));
cphi   = cos(ea_rad(3));
sphi   = sin(ea_rad(3));

T  = [ctheta * cpsi, sphi * stheta * cpsi - spsi * cphi, cphi * cpsi * stheta + sphi * spsi; ...
    ctheta * spsi, sphi * stheta * spsi + cpsi * cphi, cphi * spsi * stheta - sphi * cpsi; ...
    -stheta,                      sphi * ctheta,                      cphi * ctheta];

end


function [x, valid] = sin_cos_solve(lhs_s, lhs_c, rhs)
%
% Solves "lhs_s * sin(x) + lhs_c * cos(x) = rhs"
% for x (2 solutions). If the equation does not
% have a solution, or has complex roots, the "valid"
% output flag will be false. We do this to allow for
% code generation with only finite & real numbers.
%

tol   = 1e3 * eps;
valid = false;
x     = [nan; nan]; % (invalid solution)

if abs(lhs_s) > tol && abs(lhs_c) > tol
    % lhs_s * sin(x) + lhs_c * cos(x) = d * sin(x + p) = rhs
    % d = sqrt(lhs_s^2 + lhs_c^2)
    % p = atan2(lhs_c, lhs_s)
    arg = rhs / sqrt(lhs_s * lhs_s + lhs_c * lhs_c);
    
    if abs(arg) <= 1.
        x1 = asin(arg);
        p  = atan2(lhs_c, lhs_s);
        
        valid = true;
        x     = [wrap_to_pi(x1 - p); wrap_to_pi(pi - x1 - p)];
        
    end
    
elseif abs(lhs_s) < tol && abs(lhs_c) > tol
    % b * cos(x) = c
    arg = rhs / lhs_c;
    if abs(arg) <= 1.
        x1    = acos(arg);
        
        valid = true;
        x     = [wrap_to_pi(x1); wrap_to_pi(-x1)];
        
    end
    
elseif abs(lhs_s) > tol && abs(lhs_c) < tol
    % lhs_s * sin(x) = rhs
    arg = rhs / lhs_s;
    if abs(arg) <= 1.
        x1    = asin(arg);
        
        valid = true;
        x     = [wrap_to_pi(x1); wrap_to_pi(pi - x1)];
        
    end
    
end

% if we've arrived here with valid = false, the equation
% was either "0 = rhs" or had complex roots... invalid
% solution

end


function x = wrap_to_pi(x)
%
% Wrap angle to [-pi, pi).
%

x = mod(x + pi, 2 * pi) - pi;

end


function x = wrap_to_180(x)
%
% Wrap angle to [-180, 180).
%

x = mod(x + 180, 360) - 180;

end


function S = try_rmfield(S, fieldname)
%
% Attempts to remove a field from a
% structure.
%

if isfield(S, fieldname)
    S = rmfield(S, fieldname);
end

end


function S = merge_structs(S, S_other)
%
% Merges two structs that MUST have distinct fields.
%

% check that the fieldnames are different
fnames_other = fieldnames(S_other)';
assert(~any(isfield(S, fnames_other)));

% merge the structs
for f = fnames_other
    S.(f{:}) = S_other.(f{:});
end

end