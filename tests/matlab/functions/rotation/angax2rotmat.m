function T = angax2rotmat(angle_rad, axis)
%
% (angle [rad], axis) pair to rotation matrix. The output
% rotation matrix transforms from the rotated axes to the
% original ones: x_original = T * x_rotated, note that the
% input axis remains invariant.
%

% validate inputs
assert(isscalar(angle_rad))
assert(numel(axis) == 3)


% scalar version of the algorithm
axis = axis / max(eps, norm(axis));

c = cos(angle_rad);
s = sin(angle_rad);
C = 1. - c;

xyC = axis(1) * axis(2) * C;
xzC = axis(1) * axis(3) * C;
yzC = axis(2) * axis(3) * C;

xs = axis(1) * s;
ys = axis(2) * s;
zs = axis(3) * s;

T = [axis(1) * axis(1) * C + c, xyC - zs, xzC + ys; ...
    xyC + zs, axis(2) * axis(2) * C + c, yzC - xs; ...
    xzC - ys, yzC + xs, axis(3) * axis(3) * C + c];


% Euler-Rodrigues formula
%cAx = crossmatrix(axis / max(eps, norm(axis)));
%T   = eye(3) + sin(angle_rad) * cAx + (1 - cos(angle_rad)) * cAx * cAx;

end


%function S = crossmatrix(s)
%%
%% Returns the 3 x 3 skew-symmetric cross-product matrix of a
%% 3-component vector v, such that "cross(s, v) == S * v".
%%
%
%S = [0., -s(3), s(2); ...
%    s(3), 0., -s(1); ...
%    -s(2), s(1), 0.];
%
%end