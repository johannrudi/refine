%example_int_2d  Evaluation of integral of three 2-D functions.
%  TODO
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

%% Set MATLAB Environment

% clearing
clear

% initialize refine
refine_init


%% Set Parameters

% set mask coefficients
a_0 = [1 1 ; 1 1];
a_1 = [
    0.125  0.25   0.125  0.0   ; ...
    0.25   0.625  0.5    0.125 ; ...
    0.125  0.5    0.625  0.25  ; ...
    0.0    0.125  0.25   0.125 ...
];
a_2 = a_1;
mask = {a_0, a_1, a_2};

% set derivatives matrix with
% `num_rows = length(mask)` and `num_cols = dimension`
deriv = [0 0 ; 1 0 ; 1 0];

% set number of refine steps
numRefineSteps = 0;

% set first index of mask coefficients as a column vector with
% `num_rows = length(mask)` and `num_cols = dimension`
firstMaskIndex = [0 0 ; 1 2 ; 3 4];


%% Evaluate Refinable Integral

% compute refinable integral
[x,y] = refine(mask, deriv, numRefineSteps, ...
               'FirstMaskIndexStart',firstMaskIndex, 'Display','on');
