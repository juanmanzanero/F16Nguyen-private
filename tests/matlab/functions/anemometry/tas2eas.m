function EAS_mps = tas2eas(TAS_mps, geopotential_altitude_m, ...
    temperature_offset_K)
%
% EAS (equivalent airspeed) from TAS (true airspeed) and geopotential
% altitude, with optional temperature offset. All in SI units.
%

% zero altitude solution
tol = 1e3 * eps;
if abs(geopotential_altitude_m) < tol
    EAS_mps = TAS_mps;
    return

end


% validate inputs
if nargin < 3
    temperature_offset_K = 0.;
end


% do the thing
EAS_mps = TAS_mps * sqrt(density_sigma(geopotential_altitude_m, ...
    temperature_offset_K));

end