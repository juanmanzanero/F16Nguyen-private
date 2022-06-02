function dlreset_F16_Nguyen_clib(libalias, clib, header)
%
% Unloads & reloads the "F16_Nguyen_clib".
%

if libisloaded(libalias)
    unloadlibrary(libalias);
end

while ~libisloaded(libalias)
    loadlibrary(clib, header, 'alias', libalias);
end

end