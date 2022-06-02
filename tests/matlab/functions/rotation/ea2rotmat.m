function T = ea2rotmat(ea_rad)
%
% Rotation matrix from 'ZYX' (yaw-pitch-roll) Euler angle sequence [rad].
% The rotation matrix transforms from the rotated axes to the original
% ones: x_original = T * x_rotated.
%

% validate inputs
assert(numel(ea_rad) == 3);


cpsi   = cos(ea_rad(1));
spsi   = sin(ea_rad(1));
ctheta = cos(ea_rad(2));
stheta = sin(ea_rad(2));
cphi   = cos(ea_rad(3));
sphi   = sin(ea_rad(3));

T  = [ctheta * cpsi, sphi * stheta * cpsi - spsi * cphi, cphi * cpsi * stheta + sphi * spsi; ...
      ctheta * spsi, sphi * stheta * spsi + cpsi * cphi, cphi * spsi * stheta - sphi * cpsi; ...
            -stheta,                      sphi * ctheta,                      cphi * ctheta];

end
