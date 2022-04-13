function dlclose_F16_Nguyen_clib(libalias)
%
% Unloads the "F16_Nguyen_clib".
%

while libisloaded(libalias)
    unloadlibrary(libalias);
end

end