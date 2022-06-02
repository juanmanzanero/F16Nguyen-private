function CAS_mps = tas2cas(TAS_mps, geopotential_altitude_m, ...
    temperature_offset_K)
%
% CAS (calibrated airspeed) from TAS (true airspeed) and geopotential
% altitude, with optional temperature offset. All in SI units.
%

% zero altitude solution
tol = 1e3 * eps;
if abs(geopotential_altitude_m) < tol
    CAS_mps = TAS_mps;
    return

end


% validate inputs
if nargin < 3 || isempty(temperature_offset_K)
    temperature_offset_K = 0.;
end


% do the thing
if abs(temperature_offset_K) < tol
    [~, soundspeed, ~, ~, pressure] = isa_atmos(geopotential_altitude_m);
    
    CAS_mps = cas(TAS_mps / soundspeed, pressure);

else
    [density_SL, ~, ~, ~, ...
        pressure_SL] = isa_atmos(0., temperature_offset_K);

    [~, soundspeed, ~, ~, ...
        pressure] = isa_atmos(geopotential_altitude_m, temperature_offset_K);

    CAS_mps = cas(TAS_mps / soundspeed, pressure, density_SL, pressure_SL);

end

end