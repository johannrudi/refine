%example_int_1d  Evaluation of integral of four 1-D functions.
%  TODO
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-06-29
%  ----------------------------------------------------------------------------

%% Set MATLAB Environment

% clearing
clear

% initialize refine
refine_init


%% Set Parameters

% set mask coefficients
a_0 = [1 1];
a_1 = [0.25 0.75 0.75 0.25];
a_2 = [0.25 0.75 0.75 0.25];
a_3 = [0.25 0.75 0.75 0.25];
mask = {a_0, a_1, a_2, a_3};

% set derivatives matrix with
% `num_rows = length(mask)` and `num_cols = dimension`
deriv = [0 ; 0 ; 1 ; 1];

% set number of refine steps
numRefineSteps = 0;

% set first index of mask coefficients as a column vector with
% `num_rows = length(mask)`
firstMaskIndex = [0 ; 1 ; 2 ; 3];


%% Evaluate Refinable Integral

% compute refinable integral
[x,y] = refine(mask, deriv, numRefineSteps, ...
               'FirstMaskIndex',firstMaskIndex, 'Display','on');
