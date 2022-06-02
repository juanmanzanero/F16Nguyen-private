function TAS_mps = eas2tas(EAS_mps, geopotential_altitude_m, ...
    temperature_offset_K)
%
% TAS (true airspeed) from EAS (equivalent airspeed) and geopotential
% altitude, with optional temperature offset. All in SI units.
%

% zero altitude solution
tol = 1e3 * eps;
if abs(geopotential_altitude_m) < tol
    TAS_mps = EAS_mps;
    return

end


% validate inputs
if nargin < 3
    temperature_offset_K = 0.;
end


% do the thing
TAS_mps = EAS_mps / sqrt(density_sigma(geopotential_altitude_m, ...
    temperature_offset_K));

end
