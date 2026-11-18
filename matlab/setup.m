function setup()
%SETUP puts the sim's matlab code and the built bridge on the path
%   run once per matlab session (or add it to startup.m)
%
%   the bridge has to be built first, from the repo root:
%       cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBSIM_MATLAB=ON
%       cmake --build build-release

here = fileparts(mfilename('fullpath'));
bin = fullfile(here, 'bin');

addpath(here);
addpath(bin);

if exist('bsim_mex', 'file') ~= 3
    warning('bsim:setup', ['bsim_mex isn''t built yet, so the +bsim functions won''t work.\n' ...
        'from the repo root run:\n' ...
        '  cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DBSIM_MATLAB=ON\n' ...
        '  cmake --build build-release']);
end
end
