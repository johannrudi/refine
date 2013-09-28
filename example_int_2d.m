%EXAMPLE_INT_2D Evaluation of refinable integral of three 2-D functions.
% Evaluates a refinable integral with three factors of refinable 2-D functions
% in the integral:
%
%   I(z) = \int f_0(x) * [D^(1,0) f_1](x - z_1)
%                      * [D^(1,0) f_2](x - z_2) \dx,
%
% where `f_0, f_1, f_3` are box splines and `z_1, z_2` are 2-dimensional
% integer vectors.
%
% See also: EXAMPLE_INT_1D, EXAMPLE_INT_3D, REFINE.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

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

% combine the mask coefficients into a cell-array
mask = {a_0, a_1, a_2};

% set derivatives matrix (number of rows = number of functions,
% number of columns = dimension of the functions' domain)
deriv = [0 0 ; 1 0 ; 1 0];

% set number of refine steps
numRefineSteps = 0;

% set first index for mask coefficients (same dimensions as derivative matrix)
firstMaskIndex = [0 0 ; 1 2 ; 3 4];


%% Compute Integral Evaluation

% evaluate refinable integral
[x,y] = refine(mask, deriv, numRefineSteps, ...
               'FirstMaskIndex',firstMaskIndex, 'Display','on');

