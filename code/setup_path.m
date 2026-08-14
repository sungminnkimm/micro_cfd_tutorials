%SETUP_PATH  Add the course code folders to the MATLAB path.
%
%   Run once per MATLAB session:
%       cd D:\Projects\CFD\code
%       setup_path
%
%   After this, every script and function in the course can be called by name
%   from anywhere.

here = fileparts(mfilename('fullpath'));

addpath(fullfile(here,'lib'));
d = dir(fullfile(here,'ch*'));
for k = 1:numel(d)
    if d(k).isdir
        addpath(fullfile(here,d(k).name));
    end
end

fprintf('Microscale CFD course: path configured (%d chapter folders + lib).\n', ...
        sum([d.isdir]));
clear here d k
