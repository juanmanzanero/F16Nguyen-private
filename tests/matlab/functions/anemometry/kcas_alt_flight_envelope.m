function quad_envelope = kcas_alt_flight_envelope(KCAS_limits, ZP_ft_limits, ...
    Mach_limits, temperature_offset_K)
%
% Generates a flight envelope in [KCAS, ZP_ft], but taking into account
% the cuts with the Mach limits (use [-inf, inf] to specify the absence of
% limits). Returns a "Quad_mapping" object representing the resulting
% quadrilateral region, that has local variables "chi" & "eta" \in [-1, 1].
% KCAS and ZP can be obtained from them with
% "[KCAS, ZP_ft] = quad_envelope.evaluate(chi, eta)".
%

% validate inputs
assert(numel(KCAS_limits) == 2)
assert(KCAS_limits(1) < KCAS_limits(2))
KCAS_limits = KCAS_limits(:)';

assert(numel(ZP_ft_limits) == 2)
assert(ZP_ft_limits(1) < ZP_ft_limits(2))
ZP_ft_limits = ZP_ft_limits(:)';

if nargin < 3 || isempty(Mach_limits)
    Mach_limits = [-inf, inf];

else
    assert(numel(Mach_limits) == 2)
    assert(Mach_limits(1) < Mach_limits(2))
    Mach_limits = Mach_limits(:)';

end

if nargin < 4 || isempty(temperature_offset_K)
    temperature_offset_K               = 0.;
    [density_SL, ~, ~, ~, pressure_SL] = isa_atmos();

else
    [density_SL, ~, ~, ~, pressure_SL] = ...
        isa_atmos(0., temperature_offset_K);
    
end


% bottom, right, top, left segments of the cube
bottom = [KCAS_limits(:)'; ZP_ft_limits(1), ZP_ft_limits(1)];
right  = [KCAS_limits(2), KCAS_limits(2); ZP_ft_limits(:)'];
top    = [KCAS_limits(:)'; ZP_ft_limits(2), ZP_ft_limits(2)];
left   = [KCAS_limits(1), KCAS_limits(1); ZP_ft_limits(:)'];


% calculate the cuts with the lower & upper Mach lines
% we're flight mechanics, so we know how the Mach goes in the (KCAS, ZP_ft)
% plane... a Mach can either cut left-bottom, top-bottom, left-right or
% top-right
nMach = 200;

[line_lower_Mach, cuts_lower_Mach] = Mach_curve(Mach_limits(1), nMach, ...
    bottom, right, top, left, ...
    ZP_ft_limits, ...
    density_SL, pressure_SL, temperature_offset_K);

[line_upper_Mach, cuts_upper_Mach] = Mach_curve(Mach_limits(2), nMach, ...
    bottom, right, top, left, ...
    ZP_ft_limits, ...
    density_SL, pressure_SL, temperature_offset_K);


% if there were no intersections, we have to see that the cube is contained
% within the limits. If not, the envelope is infeasible.
if isempty(cuts_lower_Mach) && isempty(cuts_upper_Mach)
    if interval_contains([line_lower_Mach(1, end), line_upper_Mach(1, end)], ...
            KCAS_limits) && ...
            interval_contains([line_lower_Mach(1, 1), line_upper_Mach(1, 1)], ...
            KCAS_limits)

        quad_envelope = Quad_mapping(bottom, right, top, left);
        return
        
    else
        error('kcas_alt_flight_envelope: infeasible envelope.');
        
    end
    
end


% generate the region
n     = nMach + 1;
left  = nan(2, n);
right = nan(2, n);
for i = 1:n
    KCAS_left  = max([KCAS_limits(1), line_lower_Mach(1, i)]);
    KCAS_right = min([KCAS_limits(2), line_upper_Mach(1, i)]);
    if KCAS_left <= KCAS_right
        left(:, i)  = [KCAS_left; line_lower_Mach(2, i)];
        right(:, i) = [KCAS_right; line_upper_Mach(2, i)];
    end
    
end

% orient them adequately, according to the criteria defined in
% file "Quad_mapping.m"
l     = ~isnan(left(1, :));
left  = fliplr(left(:, l));
r     = ~isnan(right(1, :));
right = fliplr(right(:, r));


% handle two special cases (not strictly necessary, but yields
% beter mappings):

% a) if both Mach lines cut the KCAS_limits(2) line
if ~isempty(cuts_lower_Mach)
    l = abs(cuts_lower_Mach(1, :) - KCAS_limits(2)) < eps;
else
    l = false;
end
lowMach_cuts_maxKCAS = any(l);

if ~isempty(cuts_upper_Mach)
    h = abs(cuts_upper_Mach(1, :) - KCAS_limits(2)) < eps;
else
    h = false;
end
highMach_cuts_maxKCAS = any(h);

if lowMach_cuts_maxKCAS && highMach_cuts_maxKCAS
    lowcut = cuts_lower_Mach(:, l);
    uppcut = cuts_upper_Mach(:, h);
    
    left(:, 1) = lowcut;
    bottom     = [lowcut, uppcut];
    
    right_over_uppcut = right(2, :) > bottom(2, 2);
    right             = [bottom(:, 2), right(:, right_over_uppcut)];
    
end

% b) if both Mach lines cut the KCAS_limits(1) line
if ~isempty(cuts_lower_Mach)
    l = abs(cuts_lower_Mach(1, :) - KCAS_limits(1)) < eps;
else
    l = false;
end
lowMach_cuts_minKCAS = any(l);

if ~isempty(cuts_upper_Mach)
    h = abs(cuts_upper_Mach(1, :) - KCAS_limits(1)) < eps;
else
    h = false;
end
highMach_cuts_minKCAS = any(h);

if lowMach_cuts_minKCAS && highMach_cuts_minKCAS
    lowcut = cuts_lower_Mach(:, l);
    uppcut = cuts_upper_Mach(:, h);
    
    right(:, end) = uppcut;
    top           = [lowcut, uppcut];
    
    left_under_lowcut = left(2, :) < top(2, 1);
    left              = [left(:, left_under_lowcut), top(:, 2)];

end

% standard case
if ~(lowMach_cuts_minKCAS && highMach_cuts_minKCAS) && ...
        ~(lowMach_cuts_maxKCAS && highMach_cuts_maxKCAS)
    bottom = [left(:, 1), right(:, 1)];
    top    = [left(:, end), right(:, end)];

end


% generate the quadrilateral mapping
quad_envelope = Quad_mapping(bottom, right, top, left);

end


%
% Helper functions.
%

function [line_Mach, cuts_Mach] = Mach_curve(Mach_value, nMach, ...
    bottom, right, top, left, ...
    ZP_ft_limits, ...
    density_SL, pressure_SL, temperature_offset_K)
%
% Cuts of line "Mach = Mach_value" with the cube that has [bottom, right, top,
% left] edges, whose ZP_ft is \in ZP_ft_limits and whose KCAS_limits is \in
% KCAS_limits, all occuring in the (KCAS, ZP_ft) space. We generate the
% "Mach = Mach_value" line by taking
% "nMach" number of points.
%

if ~isfinite(Mach_value)
    line_Mach = [Mach_value, Mach_value];
    cuts_Mach = [];
    return
    
end


% helper constants
mps2kn = 1. / 0.514444;
ft2m   = 0.3048;


% (KCAS, ZP_ft) coordinates of the Mach line, in the increasing KCAS
% direction (we know that they'll go in the decreasing ZP_ft
% direction)
dZ         = (ZP_ft_limits(2) - ZP_ft_limits(1)) / nMach;
ZP_ft_Mach = zeros(1, nMach + 1);
KCAS_Mach  = zeros(1, nMach + 1);
for i = 1:nMach + 1
    ZP_ft_Mach(i) = ZP_ft_limits(1) + (nMach - i + 1) * dZ;
    
    [~, ~, ~, ~, pressure] = isa_atmos(ZP_ft_Mach(i) * ft2m, temperature_offset_K);
    
    KCAS_Mach(i) = cas(Mach_value, pressure, density_SL, pressure_SL) * mps2kn;
    
end


% a Mach line can cut either 2 times or zero times, not more, not less (the
% case where it identically touches a corner is taken as zero)
[x_bottom, y_bottom, ~] = intersections(KCAS_Mach, ZP_ft_Mach, bottom(1, :), bottom(2, :), true);
[x_right, y_right, ~]   = intersections(KCAS_Mach, ZP_ft_Mach, right(1, :), right(2, :), true);
[x_top, y_top, ~]       = intersections(KCAS_Mach, ZP_ft_Mach, top(1, :), top(2, :), true);
[x_left, y_left, ~]     = intersections(KCAS_Mach, ZP_ft_Mach, left(1, :), left(2, :), true);

nb = numel(x_bottom);
if nb > 1
    error(['kcas_alt_flight_envelope::Mach_curve: something is very wrong', ...
        ' with this algorithm...'])

elseif nb == 1
    %     cutpos = floor(i_bottom);
    %     if abs(i_bottom - cutpos) > eps
    %         KCAS_Mach  = [KCAS_Mach(1:cutpos), x_bottom, KCAS_Mach(cutpos + 1:end)];
    %         ZP_ft_Mach = [ZP_ft_Mach(1:cutpos), y_bottom, ZP_ft_Mach(cutpos + 1:end)];
    %     end

end

nr = numel(x_right);
if nr > 1
    error(['kcas_alt_flight_envelope::Mach_curve: something is very wrong', ...
        ' with this algorithm...'])

elseif nr == 1
    %     cutpos = floor(i_right);
    %     if abs(i_right - cutpos) > eps
    %         KCAS_Mach  = [KCAS_Mach(1:cutpos), x_right, KCAS_Mach(cutpos + 1:end)];
    %         ZP_ft_Mach = [ZP_ft_Mach(1:cutpos), y_right, ZP_ft_Mach(cutpos + 1:end)];
    %     end

end

nt = numel(x_top);
if nt > 1
    error(['kcas_alt_flight_envelope::Mach_curve: something is very wrong', ...
        ' with this algorithm...'])

elseif nt == 1
    %     cutpos = floor(i_top);
    %     if abs(i_top - cutpos) > eps
    %         KCAS_Mach  = [KCAS_Mach(1:cutpos), x_top, KCAS_Mach(cutpos + 1:end)];
    %         ZP_ft_Mach = [ZP_ft_Mach(1:cutpos), y_top, ZP_ft_Mach(cutpos + 1:end)];
    %     end

end

nl = numel(x_left);
if nl > 1
    error(['kcas_alt_flight_envelope::Mach_curve: something is very wrong', ...
        ' with this algorithm...'])

elseif nl == 1
    %     cutpos = floor(i_left);
    %     if abs(i_left - cutpos) > eps
    %         KCAS_Mach  = [KCAS_Mach(1:cutpos), x_left, KCAS_Mach(cutpos + 1:end)];
    %         ZP_ft_Mach = [ZP_ft_Mach(1:cutpos), y_left, ZP_ft_Mach(cutpos + 1:end)];
    %     end

end

if (nb + nr + nt + nl) == 2
    % two cuts
    cuts_Mach = [x_bottom', x_right', x_top', x_left'; ...
        y_bottom', y_right', y_top', y_left'];
    
elseif (nb + nr + nt + nl) < 2
    % one cut (with a corner of the square), or zero cuts.
    cuts_Mach = [];
    
else
    % more than two cuts... impossible (?!)
    error(['kcas_alt_flight_envelope::Mach_curve: something is very wrong', ...
        ' with this algorithm...'])

end


line_Mach = [KCAS_Mach;...
    ZP_ft_Mach];

end


function ret = intervals_overlap(interval_a, interval_b)
%
% Returns true if the two input intervals overlap, false
% otherwise.
%

assert(interval_a(1) <= interval_a(2))
assert(interval_b(1) <= interval_b(2))

ret = (interval_a(1) <= interval_b(2)) && ...
    (interval_b(1) <= interval_a(2));

end


function ret = interval_contains(container_interval, contained_interval)
%
% Returns true if the first interval contains the second.
%

assert(container_interval(1) <= container_interval(2))
assert(contained_interval(1) <= contained_interval(2))

ret = contained_interval(2) <= container_interval(2) && ...
    contained_interval(1) >= container_interval(1);

end