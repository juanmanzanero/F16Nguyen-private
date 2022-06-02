%
% Test extracting the plant's properties.
%

addpath(genpath('./functions'))


% open lib
libalias = 'libF16_Nguyen_clib';

dtor = onCleanup(@()(dlclose_F16_Nguyen_clib(libalias)));

dlreset_F16_Nguyen_clib(libalias, ...
    ['../../lib/Release/', libalias], ...
    '../../include/F16_Nguyen/F16_Nguyen_clib.h');


% characterize lib & plant
libfunctions(libalias, '-full')

plant_properties = F16_Nguyen_plant_properties(libalias);


% finalize the lib
clear dtor