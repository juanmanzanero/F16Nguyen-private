function [trim_success, trim_outputs] = F16_Nguyen_trim_plant(trim_inflags, ...
    trim_inputs, ...
    libalias, cplant, plant_properties)
%
% Evaluates the "F16_Nguyen_clib_trim_plant" function,
% at the given "x0" (linearization states), and "u0"
% (linearization inputs).
%

% validate inputs
if isstruct(trim_inflags)
    trim_inflags = cell2mat(struct2cell(trim_inflags));
    trim_inflags = int32(trim_inflags);

end

if isstruct(trim_inputs)
    trim_inputs = cell2mat(struct2cell(trim_inputs));
end


% do the thing
[err, trim_success, trim_outputs] = ...
    calllib(libalias, 'F16_Nguyen_clib_trim_plant', ...
    int32(0), ...
    zeros(plant_properties.num_trim_outputs, 1), plant_properties.num_trim_outputs, ...
    cplant, ...
    trim_inflags, plant_properties.num_trim_inflags, ...
    trim_inputs, plant_properties.num_trim_inputs);

if err ~= 0
    error(['F16_Nguyen_trim_plant:', ...
        ' "F16_Nguyen_clib_trim_plant" returned an error (%d).'], err);

end


% convert the outputs
trim_success = logical(trim_success);
trim_outputs = cell2struct(num2cell(trim_outputs), plant_properties.trim_output_names);

end