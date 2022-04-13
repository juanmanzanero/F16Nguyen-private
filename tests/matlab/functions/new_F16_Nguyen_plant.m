function cplant = new_F16_Nguyen_plant(libalias, ...
    path_aerodataset, path_enginedataset)
%
% Creates & returns an opaque handle to an "F16_Nguyen_plant"
% instance, by calling the "F16_Nguyen_clib".
%

[err, cplant] = calllib(libalias, 'F16_Nguyen_clib_new_plant', libpointer, ...
    path_aerodataset, path_enginedataset);

if err ~= 0
    error(['create_F16_Nguyen_plant:', ...
        ' "F16_Nguyen_clib_new_plant" returned an error (%d).'], err);

end

end