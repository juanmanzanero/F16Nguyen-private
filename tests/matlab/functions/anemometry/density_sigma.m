function sigma = density_sigma(geopotential_altitude_m, ...
    temperature_offset_K)
%
% Returns the (density / density_SL) ratio, given the
% geopotential altitude and an optional temperature
% offset. All in SI units.
%

% % zero altitude solution
tol = 1e3 * eps;
if abs(geopotential_altitude_m) < tol
    sigma = 1.;
    return
    
end


% validate inputs
if nargin < 2 || isempty(temperature_offset_K)
    temperature_offset_K = 0.;
end


% do the thing
if abs(temperature_offset_K) < tol
    density_SL = isa_atmos();
    density    = isa_atmos(geopotential_altitude_m);
    sigma      = density  / density_SL;
    
else
    density_SL = isa_atmos(0., temperature_offset_K);
    density    = isa_atmos(geopotential_altitude_m, temperature_offset_K);
    sigma      = density / density_SL;
    
end

end