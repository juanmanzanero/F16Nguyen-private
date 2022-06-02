function CAS_mps = eas2cas(EAS_mps, geopotential_altitude_m, ...
    temperature_offset_K)
%
% CAS (calibrated airspeed) from EAS (equivalent airspeed) and geopotential
% altitude, with optional temperature offset. All in SI units.
%

% zero altitude solution
tol = 1e3 * eps;
if abs(geopotential_altitude_m) < tol
    CAS_mps = EAS_mps;
    return

end


% validate inputs
if nargin < 3 || isempty(temperature_offset_K)
    temperature_offset_K = 0.;
end


% do the thing
if abs(temperature_offset_K) < tol
    density_SL = isa_atmos();
    
    [density, soundspeed, ~, ~, ...
        pressure] = isa_atmos(geopotential_altitude_m);
    
    CAS_mps = cas(EAS_mps / sqrt(density / density_SL) / soundspeed, ...
        pressure);
    
else
    [density_SL , ~, ~, ~, ...
        pressure_SL] = isa_atmos(0., temperature_offset_K);
    
    [density, soundspeed, ~, ~, ...
        pressure] = isa_atmos(geopotential_altitude_m, temperature_offset_K);
    
    CAS_mps = cas(EAS_mps / sqrt(density / density_SL) / soundspeed, ...
        pressure, ...
        density_SL, pressure_SL);
    
end

end