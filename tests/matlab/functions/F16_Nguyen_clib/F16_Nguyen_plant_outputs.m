function y = F16_Nguyen_plant_outputs(x, u, ...
    libalias, cplant, plant_properties)
%
% Evaluates the "F16_Nguyen_plant"'s "outputs"
% state-space function, for the given "x" (states)
% and "u" (inputs).
%

% validate inputs
if isstruct(x)
    x = cell2mat(struct2cell(x));
end

if isstruct(u)
    u = cell2mat(struct2cell(u));
end


% do the thing
[err, y] = calllib(libalias, 'F16_Nguyen_clib_plant_outputs', ...
    zeros(plant_properties.num_outputs, 1), plant_properties.num_outputs, ...
    cplant, ...
    x, plant_properties.num_states, ...
    u, plant_properties.num_inputs);

if err ~= 0
    error(['F16_Nguyen_plant_outputs:', ...
        ' "F16_Nguyen_clib_plant_outputs" returned an error (%d).'], err);

end

end