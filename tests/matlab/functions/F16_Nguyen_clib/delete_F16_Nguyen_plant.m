function delete_F16_Nguyen_plant(libalias, cplant)
%
% Deletes the opaque handle to an "F16_Nguyen_plant"
% instance (obtained by calling "new_F16_Nguyen_plant").
%

err = calllib(libalias, 'F16_Nguyen_clib_delete_plant', cplant);

if err ~= 0
    error(['delete_F16_Nguyen_plant:', ...
        ' "F16_Nguyen_clib_delete_plant" returned an error (%d).'], err);

end

end