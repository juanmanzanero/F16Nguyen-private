function KCAS = Mach2kcas(Mach, geopotential_altitude_ft, ...
    temperature_offset_K)
%
% Mach to KCAS conversion, with geopotential altitude in ft, and optional
% temperature offset. I ended up writing this one because I'm needing it
% all the time...
%

% validate inputs
if nargin < 3 || isempty(temperature_offset_K)
    temperature_offset_K = 0.;
end


% helper constants
ft2m   = 0.3048;
mps2kn = 1. / 0.514444;


% calculate the calibrated airseed in knots
ZP_m = ft2m * geopotential_altitude_ft;
tol  = 1e3 * eps;
if abs(temperature_offset_K) < tol
    [~, ~, ~, ~, pressure] = isa_atmos(ZP_m);

    KCAS = mps2kn * cas(Mach, pressure);

else
    [density_SL, ~, ~, ~, pressure_SL] = isa_atmos(0., temperature_offset_K);
    [~, ~, ~, ~, pressure]             = isa_atmos(ZP_m, temperature_offset_K);

    KCAS = mps2kn * cas(Mach, pressure, ...
        density_SL, ...
        pressure_SL);

end

end