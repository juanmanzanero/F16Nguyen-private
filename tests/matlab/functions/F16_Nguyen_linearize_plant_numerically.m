function [A, B, C, D, ...
    xdot, y, x, u] = F16_Nguyen_linearize_plant_numerically(x0, u0, ...
    deltas_numjac_states, deltas_numjac_inputs, ...
    libalias, cplant, plant_properties)
%
% Evaluates the "F16_Nguyen_clib_linearize_plant_numerically"
% function, at the given "x0" (linearization states), and "u0"
% (linearization inputs), calculating a numerical jacobian
% with central finite differences, taking "deltas_numjac_states"
% as step sizes for the states, and "deltas_numjac_inputs" as
% step sizes for the inputs.
%

% validate inputs
if isstruct(x0)
    x0 = cell2mat(struct2cell(x0));
end

if isstruct(u0)
    u0 = cell2mat(struct2cell(u0));
end

if isstruct(deltas_numjac_states)
    deltas_numjac_states = cell2mat(struct2cell(deltas_numjac_states));
end

if isstruct(deltas_numjac_inputs)
    deltas_numjac_inputs = cell2mat(struct2cell(deltas_numjac_inputs));
end


% do the thing
[err, A, B, C, D, ...
    xdot, y, x, u] = ...
    calllib(libalias, 'F16_Nguyen_clib_linearize_plant_numerically', ...
    zeros(plant_properties.num_states), plant_properties.num_states * plant_properties.num_states, ...
    zeros(plant_properties.num_states, plant_properties.num_inputs), plant_properties.num_states * plant_properties.num_inputs, ...
    zeros(plant_properties.num_outputs, plant_properties.num_states), plant_properties.num_outputs * plant_properties.num_states, ...
    zeros(plant_properties.num_outputs, plant_properties.num_inputs), plant_properties.num_outputs * plant_properties.num_inputs, ...
    zeros(plant_properties.num_states, 1), plant_properties.num_states, ...
    zeros(plant_properties.num_outputs, 1), plant_properties.num_outputs, ...
    zeros(plant_properties.num_states, 1), plant_properties.num_states, ...
    zeros(plant_properties.num_inputs, 1), plant_properties.num_inputs, ...
    cplant, ...
    x0, plant_properties.num_states, ...
    u0, plant_properties.num_inputs, ...
    deltas_numjac_states, numel(deltas_numjac_states), ...
    deltas_numjac_inputs, numel(deltas_numjac_inputs));

if err ~= 0
    error(['F16_Nguyen_linearize_plant_numerically:', ...
        ' "F16_Nguyen_clib_linearize_plant_numerically" returned an error (%d).'], err);

end

end