%refine_init  Initialize the Refine program.
%  Performs operations that are necessary in order to use the Refine program.
%  Needs to be executed once before running the program.
%
%  Syntax:
%  refine_init
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

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
