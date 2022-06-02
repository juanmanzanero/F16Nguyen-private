function kinematic_outputs = flat_Earth_kinematics(velocity_and_acceleration_are_airspeeds, ...
    acceleration_Earthaxes_mps2, velocity_Earthaxes_mps, ...
    yprdot_degps, ypr_deg, ...
    wind_acceleration_Earthaxes_mps2, wind_velocity_Earthaxes_mps)
%
% Calculates some useful kinematic variables from the input
% aircraft's acceleration, velocity and attitude, plus the
% mean wind's acceleration and velocity.
%

% helper constants
deg2rad = pi / 180.;
rad2deg = 1. / deg2rad;
g0_mps2 = 9.80665;


% validate inputs, we'll work with columns
assert(numel(acceleration_Earthaxes_mps2) == 3)
acceleration_Earthaxes_mps2 = acceleration_Earthaxes_mps2(:);

assert(numel(velocity_Earthaxes_mps) == 3)
velocity_Earthaxes_mps = velocity_Earthaxes_mps(:);

assert(numel(wind_acceleration_Earthaxes_mps2) == 3)
wind_acceleration_Earthaxes_mps2 = wind_acceleration_Earthaxes_mps2(:);

assert(numel(wind_velocity_Earthaxes_mps) == 3)
wind_velocity_Earthaxes_mps = wind_velocity_Earthaxes_mps(:);

assert(numel(yprdot_degps) == 3)
yprdot_degps = yprdot_degps(:);

assert(numel(ypr_deg) == 3)
ypr_deg = ypr_deg(:);


% input velocity and acceleration can be either
% "with-respect-to-Earth" or "with-respect-to-air"
if velocity_and_acceleration_are_airspeeds
    axyz_Earthaxes_mps2 = acceleration_Earthaxes_mps2 + wind_acceleration_Earthaxes_mps2;
    uvw_Earthaxes_mps   = velocity_Earthaxes_mps + wind_velocity_Earthaxes_mps;

else
    axyz_Earthaxes_mps2 = acceleration_Earthaxes_mps2;
    uvw_Earthaxes_mps   = velocity_Earthaxes_mps;

end


% true airspeed
VTAS_Earthaxes_mps      = uvw_Earthaxes_mps - wind_velocity_Earthaxes_mps;
dVTAS_dt_Earthaxes_mps2 = axyz_Earthaxes_mps2 - wind_acceleration_Earthaxes_mps2;
TAS_mps                 = norm(VTAS_Earthaxes_mps);


% body axes
ypr_rad      = deg2rad * ypr_deg;
T_Earth2body = ea2rotmat(ypr_rad)';
uvw_mps      = T_Earth2body * uvw_Earthaxes_mps;
pqr_radps    = angular_kinematic_relationships(deg2rad * yprdot_degps, ypr_rad);
axyz_mps2    = T_Earth2body * axyz_Earthaxes_mps2;
uvwdot_mps2  = axyz_mps2 - cross(pqr_radps, uvw_mps); % this one's called the "pseudo-acceleration"


% wind axes
VTAS_mps      = T_Earth2body * VTAS_Earthaxes_mps;
dVTAS_dt_mps2 = T_Earth2body * dVTAS_dt_Earthaxes_mps2;
VTASdot_mps2  = dVTAS_dt_mps2 - cross(pqr_radps, VTAS_mps);

[alphadot_radps, betadot_radps, ...
    alpha_rad, beta_rad, T_body2wind, ...
    pqr_bodyRwind_radps, pqr_windREarth_radps, ...
    chidot_gammadot_mudot_radps, ...
    chi_gamma_mu_rad, T_wind2Earth] = wind_axes_kinematics(VTASdot_mps2, ...
    VTAS_mps, pqr_radps, T_Earth2body);


% ground-velocity axes
[alphagrounddot_radps, betagrounddot_radps, ...
    alphaground_rad, betaground_rad, T_body2groundvel, ...
    pqr_bodyRgroundvel_radps, pqr_groundvelREarth_radps, ...
    chigrounddot_gammagrounddot_mugrounddot_radps, ...
    chiground_gammaground_muground_rad, T_groundvel2Earth] = ...
    wind_axes_kinematics(uvwdot_mps2, uvw_mps, pqr_radps, T_Earth2body);


% rate of climb & groundspeed
ROC_mps = -uvw_Earthaxes_mps(3);
GS_mps  = norm(uvw_Earthaxes_mps);


% turn radius (for the torsion we would need
% "axyzdot_Earthaxes_mps3" as an input, it seems
% too much for our purposes, but I leave it here
% nonetheless)
whatever_jerk_mps3 = [0; 0; 0];

[ground_turn_radius_m, ~] = curvature_radius_and_torsion(whatever_jerk_mps3, ...
    axyz_mps2, uvw_mps);

[air_turn_radius_m, ~] = curvature_radius_and_torsion(whatever_jerk_mps3, ...
    dVTAS_dt_mps2, VTAS_mps);

% % alternative & equivalent formulas for the turn radius
% ground_turn_radius_m = GS_mps / ...
%     sqrt(chigrounddot_gammagrounddot_mugrounddot_radps(1)^2 * cos(chiground_gammaground_muground_rad(2))^2 + ...
%     chigrounddot_gammagrounddot_mugrounddot_radps(2)^2);
% 
% air_turn_radius_m = TAS_mps / ...
%     sqrt(chidot_gammadot_mudot_radps(1)^2 * cos(chi_gamma_mu_rad(2))^2 + ...
%     chidot_gammadot_mudot_radps(2)^2);


% body-axes load factors (by the flat Earth model, g_Earthaxes = [0; 0; g0]),
% wind-axes load factors and ground-velocity-axes load factors
Nxyz_g          = T_Earth2body * ([0; 0; g0_mps2] - axyz_Earthaxes_mps2) / g0_mps2;
Nwxyz_g         = T_body2wind * Nxyz_g;
Ngroundvelxyz_g = T_body2groundvel * Nxyz_g;


% gather outputs
kinematic_outputs = struct('axyz_Earthaxes_mps2', axyz_Earthaxes_mps2, ...
    'uvw_Earthaxes_mps', uvw_Earthaxes_mps, ...
    'VTAS_Earthaxes_mps', VTAS_Earthaxes_mps, ...
    'dVTAS_dt_Earthaxes_mps2', dVTAS_dt_Earthaxes_mps2, ...
    'TAS_mps', TAS_mps, ...
    ...
    'uvw_mps', uvw_mps, ...
    'pqr_degps', rad2deg * pqr_radps, ...
    'uvwdot_mps2', uvwdot_mps2, ...
    'axyz_mps2', axyz_mps2, ...
    ...
    'VTAS_mps', VTAS_mps, ...
    'dVTAS_dt_mps2', dVTAS_dt_mps2, ...
    'VTASdot_mps2', VTASdot_mps2, ...
    'alpha_deg', rad2deg * alpha_rad, ...
    'beta_deg', rad2deg * beta_rad, ...
    'T_body2wind', T_body2wind, ...
    'alphadot_degps', rad2deg * alphadot_radps, ...
    'betadot_degps', rad2deg * betadot_radps, ...
    'pqr_bodyRwind_degps', rad2deg * pqr_bodyRwind_radps, ...
    'pqr_windREarth_degps', rad2deg * pqr_windREarth_radps, ...
    'chi_gamma_mu_deg', rad2deg * chi_gamma_mu_rad, ...
    'T_wind2Earth', T_wind2Earth, ...
    'chidot_gammadot_mudot_degps', rad2deg * chidot_gammadot_mudot_radps, ...
    ...
    'alphaground_deg', rad2deg * alphaground_rad, ...
    'betaground_deg', rad2deg * betaground_rad, ...
    'T_body2groundvel', T_body2groundvel, ...
    'alphagrounddot_degps', rad2deg * alphagrounddot_radps, ...
    'betagrounddot_degps', rad2deg * betagrounddot_radps, ...
    'pqr_bodyRgroundvel_degps', rad2deg * pqr_bodyRgroundvel_radps, ...
    'pqr_groundvelREarth_degps', rad2deg * pqr_groundvelREarth_radps, ...
    'chiground_gammaground_muground_deg', rad2deg * chiground_gammaground_muground_rad, ...
    'T_groundvel2Earth', T_groundvel2Earth, ...
    'chigrounddot_gammagrounddot_mugrounddot_degps', rad2deg * chigrounddot_gammagrounddot_mugrounddot_radps, ...
    ...
    'ROC_mps', ROC_mps, ...
    'GS_mps', GS_mps, ...
    ...
    'ground_turn_radius_m', ground_turn_radius_m, ...
    'air_turn_radius_m', air_turn_radius_m, ...
    ...
    'Nxyz_g', Nxyz_g, ...
    'Nwxyz_g', Nwxyz_g, ...
    'Ngroundvelxyz_g', Ngroundvelxyz_g);

end


%
% Helper functions.
%

function [alphadot_radps, betadot_radps, ...
    alpha_rad, beta_rad, T_body2wind, ...
    pqr_bodyRwind_radps, pqr_windREarth_radps, ...
    chidot_gammadot_mudot_radps, ...
    chi_gamma_mu_rad, T_wind2Earth] = wind_axes_kinematics(VTASdot_mps2, ...
    VTAS_mps, pqr_radps, T_Earth2body)
%
% Returns useful kinematic variables that concern the wind-axes frame.
%

%  body-to-wind rot. matrix, aoa, aos, and their
% time derivatives
[alpha_rad, beta_rad, ...
    T_body2wind] = vtas2wa(VTAS_mps);

[alphadot_radps, betadot_radps] = ...
    wind_angles_time_derivatives(VTASdot_mps2, VTAS_mps);


% wind-axes angular velocities
[pqr_bodyRwind_radps, pqr_windREarth_radps] = ...
    wind_axes_angular_velocities(pqr_radps, ...
    alphadot_radps, betadot_radps, T_body2wind);


% wind-to-Earth rot. matrix, flight-path-airspeed angles
% and their time derivatives
T_wind2Earth     = (T_body2wind * T_Earth2body)';
chi_gamma_mu_rad = rotmat2ea(T_wind2Earth);

chidot_gammadot_mudot_radps = ...
    inverse_angular_kinematic_relationships(pqr_windREarth_radps, ...
    chi_gamma_mu_rad);

end


function [alpha_rad, beta_rad, T_body2wind] = vtas2wa(VTAS_mps)
%
% Calculates the body axes' "wind angles" of (attack, sideslip) [rad]
% from input "VTAS_mps" (true airspeed vector, projected onto these body
% axes and in [m / s]). Also returns "T_body2wind", i.e., the rot. matrix
% from body to wind axes.
%

% if norm(VTAS_mps) > eps
    % get the speed components
    w = VTAS_mps(3);
    v = VTAS_mps(2);
    u = VTAS_mps(1);
    
    % use a Fortran-style sign (Matlab's sign returns 0 for sign(0))
    sgn_u = sign(u) + (u == 0);
    
    % calculate alpha \in [-pi / 2, pi / 2]
    salpha    = sgn_u * w;
    calpha    = abs(u); % always positive
    alpha_rad = atan2(salpha, calpha);

    % calculate beta \in [-pi, pi)
    sbeta    = v;
    cbeta    = sgn_u * sqrt(u * u + w * w);
    beta_rad = mod(atan2(sbeta, cbeta) + pi, 2 * pi) - pi; % atan2 returns an angle in [-pi, pi], this takes it to [-pi, pi)

    % calculate the rot. matrix
    calpha      = cos(alpha_rad);
    salpha      = sin(alpha_rad);
    cbeta       = cos(beta_rad);
    sbeta       = sin(beta_rad);
    T_body2wind = [calpha * cbeta, sbeta,  salpha * cbeta; ...
        -calpha * sbeta, cbeta, -salpha * sbeta; ...
        -salpha,     0,          calpha];

% else
%     % agreed static solution
%     alpha_rad   = 0.;
%     beta_rad    = 0.;
%     T_body2wind = eye(3);
%     
% end

end


function [alphadot_radps, betadot_radps] = ...
    wind_angles_time_derivatives(VTASdot_mps2, VTAS_mps)
%
% Calculates the time derivatives of the body axes' "wind angles" of
% (attack, sideslip) [rad / s] from inputs "VTASdot_mps2" (airspeed
% pseudo-acceleration) and "VTAS_mps" (true airspeed vector), both
% projected onto these body axes.
%

% get the speed & pseudoacceleration components
w = VTAS_mps(3);
v = VTAS_mps(2);
u = VTAS_mps(1);

wdot = VTASdot_mps2(3);
vdot = VTASdot_mps2(2);
udot = VTASdot_mps2(1);


% use a Fortran-style sign (Matlab's sign returns 0 for sign(0))
sgn_u = sign(u) + (u == 0);


% calculate some helper airspeed norms and the desired time
% derivatives
TAS = norm(VTAS_mps);
vl2 = u * u + w * w;

alphadot_radps = (wdot * u - udot * w) / max(vl2, eps);

betadot_radps = sgn_u * (vdot * vl2 - v * (udot * u + wdot * w)) / ...
    max(TAS * TAS * sqrt(vl2), eps);

end


function [pqr_bodyRwind_radps, pqr_windREarth_radps] = ...
    wind_axes_angular_velocities(pqr_radps, ...
    alphadot_radps, betadot_radps, T_body2wind)
%
% Returns two angular velocity vectors, in [rad / s]:
% "pqr_windREarth_radps", i.e., the angular velocity
% of the wind axes relative to the Earth axes frame,
% projected on wind axes; and "pqr_bodyRwind_radps",
% i.e., the angular velocity of the body axes relative
% to the wind axes frame, projected on body axes.
%

% this one is just "-betadot * z_wind + alphadot * y_body"
salpha              = -T_body2wind(3, 1);
calpha              = T_body2wind(3, 3);
pqr_bodyRwind_radps = [betadot_radps * salpha; ...
    alphadot_radps; ...
    -betadot_radps * calpha];


% this one has the famous [pw, qw, rw] components
pqr_windREarth_radps = T_body2wind * ...
    (pqr_radps - pqr_bodyRwind_radps);

end


function pqr_radps = angular_kinematic_relationships(yprdot_radps, ypr_rad)
%
% Returns the angular velocity of a body [rad / s], projected onto its
% body axes frame, calculated from its Euler angles (in 'ZYX' ->
% yaw-pitch-roll sequence) and their time derivatives.
%

psidot_ctheta = yprdot_radps(1) * cos(ypr_rad(2));
cphi          = cos(ypr_rad(3));
sphi          = sin(ypr_rad(3));

pqr_radps = [yprdot_radps(3) - yprdot_radps(1) * sin(ypr_rad(2)); ...
    yprdot_radps(2) * cphi + psidot_ctheta * sphi; ...
    -yprdot_radps(2) * sphi + psidot_ctheta * cphi];

end


function yprdot_radps = ...
    inverse_angular_kinematic_relationships(pqr_radps, ypr_rad)
%
% Returns the time derivatives of the Euler angles of a body
% (in 'ZYX' -> yaw-pitch-roll sequence) [rad / s], calculated
% from its body axes angular velocity and its attitude.
%

cphi      = cos(ypr_rad(3));
sphi      = sin(ypr_rad(3));
omega_lat = pqr_radps(2) * sphi + pqr_radps(3) * cphi;

yprdot_radps = [omega_lat / cos(ypr_rad(2)); ...
    pqr_radps(2) * cphi - pqr_radps(3) * sphi; ...
    pqr_radps(1) + omega_lat * tan(ypr_rad(2))];

end


function [curvature_radius_m, torsion_1pm] = ...
    curvature_radius_and_torsion(jerk_mps3, ...
    acceleration_mps2, velocity_mps)
%
% Returns the curvature radius [m] and torsion [1 / m] of a
% trajectory, from velocity, acceleration and jerk (i.e.,
% the time derivative of the acceleration). Note that
% "curvature = 1. / curvature_radius".
%

normcross_v_a = norm(cross(velocity_mps, ...
    acceleration_mps2));

curvature_radius_m = norm(velocity_mps)^3 / normcross_v_a;

torsion_1pm = dot(velocity_mps, ...
    cross(acceleration_mps2, jerk_mps3)) / ...
    normcross_v_a^2;

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


function [ea_rad, gymbal_lock] = rotmat2ea(T, gymbal_lock_pitch_tol_deg)
%
% Euler angles in 'ZYX' sequence [rad] from rotation matrix. The rotation
% matrix transforms from the rotated axes to the original
% ones, i.e., "x_original = T * x_rotated". Returns an additional flag
% signalling gymbal lock.
%

if nargin < 2 || isempty(gymbal_lock_pitch_tol_deg)
    gymbal_lock_pitch_tol_deg = 1e-3; % deg tolerance for gymbal lock
else
    assert(isscalar(gymbal_lock_pitch_tol_deg))
end


ea_rad = zeros(3, 1);

% pitch
gymbal_lock_zone = min(1., max(0., cosd(gymbal_lock_pitch_tol_deg)));
sinp             = -T(3, 1);
if abs(sinp) <= gymbal_lock_zone
    gymbal_lock = false;
    ea_rad(2)   = asin(sinp);
    
else
    % gymbal lock, pitch = \pm (pi * 0.5)
    gymbal_lock = true;
    ea_rad(2)   = sign(sinp) * pi * 0.5;
    
end

% yaw
ea_rad(1) = atan2(T(2, 1), T(1, 1));

% roll
ea_rad(3) = atan2(T(3, 2), T(3, 3));

% wrap to [-pi, pi)
ea_rad = mod(ea_rad + pi, 2 * pi) - pi;

end