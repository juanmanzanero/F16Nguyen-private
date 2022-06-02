function EAS_mps = cas2eas(CAS_mps, geopotential_altitude_m, ...
    temperature_offset_K)
%
% EAS (equivalent airspeed) from CAS (calibrated airspeed) and geopotential
% altitude, with optional temperature offset. All in SI units.
%

% zero altitude solution
tol = 1e3 * eps;
if abs(geopotential_altitude_m) < tol
    EAS_mps = CAS_mps;
    return

end


% validate inputs
if nargin < 3
    temperature_offset_K = 0.;
end


% do the thing
tol  = 1e-8;
seed = CAS_mps;
fun  = @(eas)(eas2cas(eas, geopotential_altitude_m, temperature_offset_K) - CAS_mps);

[EAS_mps, ~, exitflag] = fzero(fun, seed, struct('TolX', tol));

if exitflag ~= 1
    warning('cas2eas: the root finding algorithm failed, result will be approximate.');
    
    EAS_mps = seed;
    
end

end