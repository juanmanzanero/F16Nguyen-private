%
% Tests a plant linearization, at a S & L trim point.
%

addpath('./functions');
addpath(genpath('/home/c83833/projects/arbs/boom3d/visualize/functions'));


% helper constants
ft2m    = 0.3048;
kn2mps  = 0.514444;
deg2rad = pi / 180.;


% open lib
libalias = 'libF16_Nguyen_clib';

dtor = onCleanup(@()(dlclose_F16_Nguyen_clib(libalias)));

dlreset_F16_Nguyen_clib(libalias, ...
    ['../../lib/Release/', libalias], ...
    '../../include/F16_Nguyen/F16_Nguyen_clib.h');


% create the plant
plant_properties = F16_Nguyen_plant_properties(libalias);

cplant = new_F16_Nguyen_plant(libalias, ...
    '../../datasets/aero_betasym', '../../datasets/engine');


% trim the plant
trim_inflags                         = plant_properties.default_trim_inflags;
trim_inflags.steady_trim             = true;
trim_inflags.straight_and_level_trim = true;
trim_inflags.disable_lef             = false;

trim_inputs                       = plant_properties.default_trim_inputs;
trim_inputs.KCAS                  = 350;
trim_inputs.ZP_ft                 = 27000;
trim_inputs.flight_path_angle_deg = 25;
trim_inputs.mass_kg               = 11e3;
trim_inputs.xcg_per_MAC           = 0.15;

[trim_success, ...
    trim_outputs] = F16_Nguyen_trim_plant(trim_inflags, trim_inputs, ...
    libalias, cplant, plant_properties);

if ~trim_success
    error('test_linearization: unsuccessful trim...')
end


% set the linearization states
x0 = cell2struct(num2cell(zeros(plant_properties.num_states, 1)), ...
    plant_properties.state_names);

ZP_m     = ft2m * trim_inputs.ZP_ft;
TAS_mps  = cas2tas(kn2mps * trim_inputs.KCAS, ZP_m);
x0.u_mps = TAS_mps * cosd(trim_outputs.aoa_deg);
x0.w_mps = TAS_mps * sind(trim_outputs.aoa_deg);

quat_body2Earth  = ea2quat(deg2rad * [0; trim_outputs.pitch_deg; 0]);
x0.qw_body2Earth = quat_body2Earth(1);
x0.qx_body2Earth = quat_body2Earth(2);
x0.qy_body2Earth = quat_body2Earth(3);
x0.qz_body2Earth = quat_body2Earth(4);

x0.zEarth_m   = -ZP_m;
x0.dh_deg     = trim_outputs.dh_deg;
x0.dlef_deg   = trim_outputs.dlef_deg;
x0.P3_percent = trim_outputs.P3_percent;


% set the linearization inputs
u0                  = plant_properties.default_inputs;
u0.dh_dmd_deg       = trim_outputs.dh_deg;
u0.dlef_dmd_deg     = trim_outputs.dlef_deg;
u0.throttle_percent = trim_outputs.throttle_percent;
u0.mass_kg          = trim_inputs.mass_kg;
u0.xcg_per_MAC      = trim_inputs.xcg_per_MAC;


% obtain the linearization
fprintf('test_linearization: autodiff linearization...')
tbegin = tic;
[A, B, C, D, ...
    xdot, y, x, u] = F16_Nguyen_linearize_plant(x0, u0, ...
    libalias, cplant, plant_properties);
fprintf(' done in %.17f\n', toc(tbegin))


% compare against a numerical linearization
delta_numjac = 1e-8;

fprintf('test_linearization: numerical linearization...')
tbegin = tic;
[A_numerical, B_numerical, C_numerical, D_numerical, ...
    xdot_numerical, y_numerical, x_numerical, u_numerical] = ...
    F16_Nguyen_linearize_plant_numerically(x0, u0, ...
    delta_numjac, delta_numjac, ...
    libalias, cplant, plant_properties);
fprintf(' done in %.17f\n', toc(tbegin))


% compare against a numerical linearization in Matlab
[A_matlab, B_matlab, C_matlab, D_matlab, ...
    xdot_matlab, y_matlab, x_matlab, u_matlab]  = ...
    linearize_plant_numerically(libalias, cplant, plant_properties, ...
    cell2mat(struct2cell(x0)), cell2mat(struct2cell(u0)), delta_numjac);


% eliminate the "*_turn_radius" outputs, since they fuck with
% our comparisons and they're inifinity anyways
turn_radius_rows                 = contains(plant_properties.output_names, 'turn_radius');
C(turn_radius_rows, :)           = [];
C_numerical(turn_radius_rows, :) = [];
C_matlab(turn_radius_rows, :)    = [];
D(turn_radius_rows, :)           = [];
D_numerical(turn_radius_rows, :) = [];
D_matlab(turn_radius_rows, :)    = [];
y(turn_radius_rows)              = [];
y_numerical(turn_radius_rows)    = [];
y_matlab(turn_radius_rows)       = [];


% compare numerical linearization results, they should be super similar...
tol = delta_numjac / 10;
assert(max(abs(A_numerical(:) - A_matlab(:))) < tol)
assert(max(abs(B_numerical(:) - B_matlab(:))) < tol)
assert(max(abs(C_numerical(:) - C_matlab(:))) < tol)
assert(max(abs(D_numerical(:) - D_matlab(:))) < tol)
assert(max(abs(xdot_numerical(:) - xdot_matlab(:))) < tol)
assert(max(abs(y_numerical(:) - y_matlab(:))) < tol)
assert(max(abs(x_numerical(:) - x_matlab(:))) < tol)
assert(max(abs(u_numerical(:) - u_matlab(:))) < tol)


% assess the autodiff vs numerical errors
[errA, rowA, colA] = matrix_relative_error(A, A_numerical);
[errB, rowB, colB] = matrix_relative_error(B, B_numerical);
[errC, rowC, colC] = matrix_relative_error(C, C_numerical);
[errD, rowD, colD] = matrix_relative_error(D, D_numerical);
[errxdot, rowxdot] = matrix_relative_error(xdot, xdot_numerical);
[erry, rowy]       = matrix_relative_error(y, y_numerical);
[errx, rowx]       = matrix_relative_error(x, x_numerical);
[erru, rowu]       = matrix_relative_error(u, u_numerical);


% summary, comparing autodiff vs numerical
fprintf('test_linearization: summary, autodiff vs numerical:\n\n')

fprintf('  -> max. rel. error A    = %.17f, at (i, j) = (%d, %d), "d(%s)/dt" vs. "%s".\n', ...
    errA, rowA, colA, plant_properties.state_names{rowA}, plant_properties.state_names{colA});

fprintf('  -> max. rel. error B    = %.17f, at (i, j) = (%d, %d), "d(%s)/dt" vs. "%s".\n', ...
    errB, rowB, colB, plant_properties.state_names{rowB}, plant_properties.input_names{colB});

fprintf('  -> max. rel. error C    = %.17f, at (i, j) = (%d, %d), "d(%s)/dt" vs. "%s".\n', ...
    errC, rowC, colC, plant_properties.output_names{rowC}, plant_properties.state_names{colC});

fprintf('  -> max. rel. error D    = %.17f, at (i, j) = (%d, %d), "d(%s)/dt" vs. "%s".\n', ...
    errD, rowD, colD, plant_properties.output_names{rowD}, plant_properties.input_names{colD});

fprintf('  -> max. rel. error xdot = %.17f, at pos    = %d, "d(%s)/dt".\n', ...
    errxdot, rowxdot, plant_properties.state_names{rowxdot});

fprintf('  -> max. rel. error y    = %.17f, at pos    = %d, "%s".\n', ...
    erry, rowy, plant_properties.output_names{rowy});
    
fprintf('  -> max. rel. error x    = %.17f, at pos    = %d, "%s".\n', ...
    errx, rowx, plant_properties.state_names{rowx});

fprintf('  -> max. rel. error u    = %.17f, at pos    = %d, "%s".\n', ...
        erru, rowu, plant_properties.input_names{rowu});


% finalize the plant and the lib
delete_F16_Nguyen_plant(libalias, cplant);
clear dtor


%
% Helper functions.
%

function [A, B, C, D, xdot, y, x, u] = ...
    linearize_plant_numerically(libalias, cplant, plant_properties, ...
    x, u, delta_numjac)
%
% Performs a numerical linearization of a plant,
% using central finite differences, at the
% linearization point defined by "[x, u]".
%

% check input sizes
assert(numel(x) == plant_properties.num_states);
assert(numel(u) == plant_properties.num_inputs);


% finite difference denominator
den = 1. / (2 * delta_numjac);


% state matrix "A = dxdot/dx"; output matrix "C = dy/dx"
A = zeros(plant_properties.num_states, plant_properties.num_states);
C = zeros(plant_properties.num_outputs, plant_properties.num_states);
for col = 1:plant_properties.num_states
    x0 = x(col);
    
    x(col)  = x0 + delta_numjac;
    y_p     = F16_Nguyen_plant_outputs(x, u, ...
        libalias, cplant, plant_properties);
    xdot_p  = F16_Nguyen_plant_derivatives(x, y_p, u, ...
        libalias, cplant, plant_properties);
    
    x(col) = x0 - delta_numjac;
    y_m    = F16_Nguyen_plant_outputs(x, u, ...
        libalias, cplant, plant_properties);
    xdot_m = F16_Nguyen_plant_derivatives(x, y_m, u, ...
        libalias, cplant, plant_properties);

    x(col) = x0;

    A(:, col) = den * (xdot_p - xdot_m);
    C(:, col) = den * (y_p - y_m);

end


% input matrix "B = dxdot/du"; feedthrough matrix "D = dy/du"
B = zeros(plant_properties.num_states, plant_properties.num_inputs);
D = zeros(plant_properties.num_outputs, plant_properties.num_inputs);
for col = 1:plant_properties.num_inputs
    u0 = u(col);
    
    u(col)  = u0 + delta_numjac;
    y_p     = F16_Nguyen_plant_outputs(x, u, ...
        libalias, cplant, plant_properties);
    xdot_p  = F16_Nguyen_plant_derivatives(x, y_p, u, ...
        libalias, cplant, plant_properties);
    
    u(col) = u0 - delta_numjac;
    y_m    = F16_Nguyen_plant_outputs(x, u, ...
        libalias, cplant, plant_properties);
    xdot_m = F16_Nguyen_plant_derivatives(x, y_m, u, ...
        libalias, cplant, plant_properties);

    u(col) = u0;

    B(:, col) = den * (xdot_p - xdot_m);
    D(:, col) = den * (y_p - y_m);

end


% outputs & statedots at the linearization point
y = F16_Nguyen_plant_outputs(x, u, ...
    libalias, cplant, plant_properties);

xdot = F16_Nguyen_plant_derivatives(x, y, u, ...
    libalias, cplant, plant_properties);

end


function [err, row, col] = matrix_relative_error(A, B)
assert(isequal(size(A), size(B)))

[err, pos] = max(abs(A(:) - B(:)));

den = min(abs(A(pos)), abs(B(pos)));
if den > 1e3 * eps
   err = err / den;
end

if numel(size(A)) == 2
    [row, col] = ind2sub(size(A), pos);

else
    row = pos;
    col = [];

end

end
