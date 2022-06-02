function [ea_rad, gymbal_lock] = rotmat2ea(T, gymbal_lock_pitch_tol_deg)
%
% Euler angles ('ZYX') [rad] from rotation matrix. The rotation matrix
% transforms from the rotated axes to the original
% ones: x_original = T * x_rotated. Returns an additional flag
% signalling gymbal lock.
%

% validate inputs
assert(all(size(T) == [3, 3]))

if nargin < 2 || isempty(gymbal_lock_pitch_tol_deg)
    gymbal_lock_pitch_tol_deg = 1e-3; % deg tolerance for gymbal lock
end


% initialize output
ea_rad = zeros(3, 1);


% pitch
gymbal_lock_zone = min(1., max(0., cosd(gymbal_lock_pitch_tol_deg)));
sinp             = -T(3, 1);
if abs(sinp) < gymbal_lock_zone
    gymbal_lock = false;
    ea_rad(2)   = asin(sinp);

else
    % gymbal lock, pitch = \pm (pi * 0.5)
    gymbal_lock = true;
    ea_rad(2)   = sign(sinp) * pi * 0.5;

end


% yaw
ea_rad(1) = atan2(T(2, 1), T(1, 1));


% roll
ea_rad(3) = atan2(T(3, 2), T(3, 3));


% wrap to [-pi, pi)
ea_rad = mod(ea_rad + pi, 2 * pi) - pi;

end
