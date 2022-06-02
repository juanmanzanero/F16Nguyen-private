function q = ea2quat(ea_rad)
%
% Quaternion from 'ZYX' (yaw-pitch-roll) Euler angle sequence [rad].
% The output quaternion transforms from the rotated axes to the
% original ones, like the input Euler angles do:
%
% [0; x_original] = quat_multiply(quat_multiply(q, [0; x_rotated]), quat_conjugate(q)).
%

% validate inputs
assert(numel(ea_rad) == 3)


midy = 0.5 * ea_rad(1);
midp = 0.5 * ea_rad(2);
midr = 0.5 * ea_rad(3);
siny = sin(midy);
sinp = sin(midp);
sinr = sin(midr);
cosy = cos(midy);
cosp = cos(midp);
cosr = cos(midr);

q = [cosy * cosp * cosr + siny * sinp * sinr; ...
     cosy * cosp * sinr - siny * sinp * cosr; ...
     cosy * sinp * cosr + siny * cosp * sinr; ...
     siny * cosp * cosr - cosy * sinp * sinr];

end