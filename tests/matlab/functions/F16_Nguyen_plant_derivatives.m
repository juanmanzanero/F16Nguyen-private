function xdot = F16_Nguyen_plant_derivatives(x, y, u, ...
    libalias, cplant, plant_properties)
%
% Evaluates the "F16_Nguyen_plant"'s "derivatives"
% state-space function, for the given "x" (states),
% "y" (outputs) and "u" (inputs).
%

% validate inputs
if isstruct(x)
    x = cell2mat(struct2cell(x));
end

if isstruct(y)
    y = cell2mat(struct2cell(y));
end

if isstruct(u)
    u = cell2mat(struct2cell(u));
end


% do the thing
[err, xdot] = calllib(libalias, 'F16_Nguyen_clib_plant_derivatives', ...
    zeros(plant_properties.num_states, 1), plant_properties.num_states, ...
    cplant, ...
    x, plant_properties.num_states, ...
    y, plant_properties.num_outputs, ...
    u, plant_properties.num_inputs);

if err ~= 0
    error(['F16_Nguyen_plant_derivatives:', ...
        ' "F16_Nguyen_clib_plant_derivatives" returned an error (%d).'], err);

end

end