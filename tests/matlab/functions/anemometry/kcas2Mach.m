function Mach = kcas2Mach(KCAS, geopotential_altitude_ft, ...
    temperature_offset_K)
%
% KCAS to Mach conversion, with geopotential altitude in ft, and optional
% temperature offset. I ended up writing this one because I'm needing it
% all the time...
%

% validate inputs
if nargin < 3
    temperature_offset_K = 0.;
end


% helper constants
ft2m   = 0.3048;
kn2mps = 0.514444;


% calculate the Mach number
ZP_m            = ft2m * geopotential_altitude_ft;
[~, soundspeed] = isa_atmos(ZP_m, temperature_offset_K);

Mach = cas2tas(kn2mps * KCAS, ZP_m, ...
    temperature_offset_K) / soundspeed;

end