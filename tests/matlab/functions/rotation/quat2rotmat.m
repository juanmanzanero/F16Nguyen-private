function T = quat2rotmat(q)
%
% Rotation matrix from quaternion. The rotation matrix transforms
% from the rotated axes to the original ones, like the quaternion
% does:
%
%    x_original      = T * x_rotated,
%    [0; x_original] = quat_multiply(quat_multiply(q, [0; x_rotated]), quat_conjugate(q)).
%
% NOTE: the input quaternion does not need to be normalized.
%

% validate inputs
assert(numel(q) == 4)


% calculate output matrix
xx  = q(2) * q(2);
yy  = q(3) * q(3);
zz  = q(4) * q(4);
two = 2. / (q(1) * q(1) + xx + yy + zz);
xx2 = xx * two;
yy2 = yy * two;
zz2 = zz * two;
wx2 = q(1) * q(2) * two;
wy2 = q(1) * q(3) * two;
wz2 = q(1) * q(4) * two;
xy2 = q(2) * q(3) * two;
xz2 = q(2) * q(4) * two;
yz2 = q(3) * q(4) * two;

T = [1. - yy2  - zz2, xy2 - wz2, xz2 + wy2; ...
    xy2 + wz2, 1. - xx2 - zz2, yz2 - wx2; ...
    xz2 - wy2, yz2 + wx2, 1. - xx2 - yy2];

end