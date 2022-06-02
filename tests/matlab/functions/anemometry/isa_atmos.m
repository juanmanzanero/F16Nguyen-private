function [rho_kgpm3, a_mps, nu_m2ps, temp_K, press_Pa] = ...
    isa_atmos(geopotential_altitude_m, temperature_offset_K)
%
% Lighter version of the 'atmos' function by  Sky Sartorius,
% "sky-s-standard-atmosphere" package, found at
% "http://www.mathworks.com/matlabcentral/fileexchange/28135"
%
% Only supports Scalar I/O, with ALL units in SI system.
%
% The ISA atmosphere model has two inputs:
%
%    * geopotential altitude [m] (= z * RE / (z + RE), z = geometric altitude, RE = 6356766 m),
%    * ISA temperature offset [K].
%
% And returns the following outputs:
%
%    * rho   = air density [kg / m^3],
%    * a     = speed of sound [m / s],
%    * nu    = kinematic viscosity [m^2 / s] ( = mu / rho),
%    * temp  = temperature [K],
%    * press = pressure [Pa].
%

% validate inputs
if nargin < 1
    % return the standard ISA Sea Level values at zero temp. offset
    % (call "isa_atmos(0., 0.)" and printf the results with '%.17g' format)
    rho_kgpm3 = 1.2249990365985557;
    a_mps     = 340.29412435568145;
    nu_m2ps   = 1.460652272024148e-05;
    temp_K    = 288.14999999999998;
    press_Pa  = 101325.;
    return

elseif nargin < 2 || isempty(temperature_offset_K)
    temperature_offset_K = 0.;

end


% helper constants
R     = 287.0531; % J / kg / K
gamma = 1.4;
g0    = 9.80665; % m / s^2
Bs    = 1.457932654517625e-06; % kg / m / s / K^(1 / 2)
S     = 110.4; % K

% lapse rate, Ki [Celsius / m]
K = [-0.0065; ... % Troposphere
    0; ... % Tropopause
    0.001; ... % Stratosphere1
    0.0028; ... % Stratosphere2
    0; ... % Stratopause
    -0.0028; ... % Mesosphere1
    -0.002; ... % Mesosphere2
    0]; % Mesopause

% base temp, Ti [K]
T = [288.15; ... % Troposphere
    216.65; ... % Tropopause
    216.65; ... % Stratosphere1
    228.65; ... % Stratosphere2
    270.65; ... % Stratopause
    270.65; ... % Mesosphere1
    214.65; ... % Mesosphere2
    186.94590831019]; % Mesopause

% base geop. alt, Hi [m]
H = [0; ... % Troposphere
    11000; ... % Tropopause
    20000; ... % Stratosphere1
    32000; ... % Stratosphere2
    47000; ... % Stratopause
    51000; ... % Mesosphere1
    71000; ... % Mesosphere2
    84852.0458449057]; % Mesopause

% base pressure, P [Pa]
P = [101325; ...  % Troposphere
    22632.0400950078; ... % Tropopause
    5474.87742428105; ... % Stratosphere1
    868.015776620216; ... % Stratosphere2
    110.90577336731; ... % Stratopause
    66.9385281211797; ... % Mesosphere1
    3.9563921603966; ... % Mesosphere2
    0.373377173762337]; % Mesopause


% calculate temperature & pressure
if geopotential_altitude_m <= H(2)
    Ton_Ti   = 1. + K(1) * (geopotential_altitude_m - H(1)) / T(1);
    temp_K   = Ton_Ti * T(1) + temperature_offset_K;
    press_Pa = P(1) * Ton_Ti^(-g0 / (K(1) * R));

elseif geopotential_altitude_m > H(2) && geopotential_altitude_m <= H(3)
    temp_K   = T(2) + temperature_offset_K;
    press_Pa = P(2) * exp(-g0 * (geopotential_altitude_m - H(2)) / (T(2) * R));
    
elseif geopotential_altitude_m > H(3) && geopotential_altitude_m <= H(4)
    Ton_Ti   = 1. + K(3) * (geopotential_altitude_m - H(3)) / T(3);
    temp_K   = Ton_Ti * T(3) + temperature_offset_K;
    press_Pa = P(3) * Ton_Ti^(-g0 / (K(3) * R));
    
elseif geopotential_altitude_m > H(4) && geopotential_altitude_m <= H(5)
    Ton_Ti   = 1. + K(4) * (geopotential_altitude_m - H(4)) / T(4);
    temp_K   = Ton_Ti * T(4) + temperature_offset_K;
    press_Pa = P(4) * Ton_Ti^(-g0 / (K(4) * R));
    
elseif geopotential_altitude_m > H(5) && geopotential_altitude_m <= H(6)
    temp_K   = T(5) + temperature_offset_K;
    press_Pa = P(5) * exp(-g0 * (geopotential_altitude_m - H(5)) / (T(5) * R));
    
elseif geopotential_altitude_m > H(6) && geopotential_altitude_m <= H(7)
    Ton_Ti   = 1. + K(6) * (geopotential_altitude_m - H(6)) / T(6);
    temp_K   = Ton_Ti * T(6) + temperature_offset_K;
    press_Pa = P(6) * Ton_Ti^(-g0 / (K(6) * R));
    
elseif geopotential_altitude_m > H(7) && geopotential_altitude_m <= H(8)
    Ton_Ti    = 1. + K(7) * (geopotential_altitude_m - H(7)) / T(7);
    temp_K   = Ton_Ti * T(7) + temperature_offset_K;
    press_Pa = P(7) * Ton_Ti^(-g0 / (K(7) * R));
    
else
    temp_K   = T(8)+ temperature_offset_K;
    press_Pa = P(8) * exp(-g0 * (geopotential_altitude_m - H(8)) / (T(8) * R));

end


% calculate the rest of variables (ideal gas + Sutherland viscosity law)
rho_kgpm3 = press_Pa / temp_K / R;
a_mps     = sqrt(gamma * R * temp_K);
nu_m2ps   = Bs * temp_K^1.5 / (temp_K + S) / rho_kgpm3;

end