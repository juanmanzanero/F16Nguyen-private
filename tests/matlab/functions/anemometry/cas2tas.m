function TAS_mps = cas2tas(CAS_mps, geopotential_altitude_m, ...
    temperature_offset_K)
%
% TAS (true airspeed) from CAS (calibrated airspeed) and geopotential
% altitude, with optional temperature offset. All in SI units.
%

% zero altitude solution
tol = 1e3 * eps;
if abs(geopotential_altitude_m) < tol
    TAS_mps = CAS_mps;
    return

end


% validate inputs
if nargin < 3
    temperature_offset_K = 0.;
end


% do the thing
tol  = 1e-8;
seed = eas2tas(CAS_mps, geopotential_altitude_m, temperature_offset_K);
fun  = @(tas)(tas2cas(tas, geopotential_altitude_m, temperature_offset_K) - CAS_mps);

[TAS_mps, ~, exitflag] = fzero(fun, seed, struct('TolX', tol));

if exitflag ~= 1
    warning('cas2tas: the root finding algorithm failed, result will be approximate.');
    
    TAS_mps = seed;
    
end

end