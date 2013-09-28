%REFINE_INIT Initialize the Refine program.
% Performs operations that are necessary in order to use the Refine program.
% Needs to be executed once per MATLAB session before running the
% main program: `refine`.
%
% Syntax:
% REFINE_INIT
%
% See also: REFINE.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function refine_init()

%% Add Directories to MATLAB's Search Path

% list of subdirectories
subdir = {'src', 'test'};

% get search path
p = path;

% add subdirectories to search path
for k = 1:length(subdir)
    s = [pwd filesep subdir{k}];
    if isempty(strfind(p, s)) % if directory is not in search path
        addpath(s);
    end
end

% end function
end
