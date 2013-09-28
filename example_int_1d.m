%EXAMPLE_INT_1D Evaluation of refinable integral of four 1-D functions.
% Evaluates a refinable integral with four factors of refinable 1-D functions
% in the integral:
%
%   I(z) = \int N_0(x) * [D^1 N_3](x - z_1)
%                      * [D^1 N_3](x - z_2)
%                      * [D^1 N_3](x - z_3) \dx,
%
% where `N_0` and `N_3` are cardinal B-Splines of order 0 and 3 (see example
% `example_fnc_1d_bspline`) and `z_1, z_2, z_3` are integers.
%
% See also: EXAMPLE_INT_2D, EXAMPLE_INT_3D, REFINE.
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

% set mask coefficients for each function in the integral
a_0 = [1 1];
a_1 = [0.25 0.75 0.75 0.25];
a_2 = [0.25 0.75 0.75 0.25];
a_3 = [0.25 0.75 0.75 0.25];

% combine the mask coefficients into a cell-array
mask = {a_0, a_1, a_2, a_3};

% set derivatives matrix (number of rows = number of functions,
% number of columns = dimension of the functions' domain)
deriv = [0 ; ...
         0 ; ...
         1 ; ...
         1];

% set number of refine steps
numRefineSteps = 0;

% set first index for mask coefficients (same dimensions as derivative matrix)
firstMaskIndex = [0 ; 1 ; 2 ; 3];


%% Compute Integral Evaluation

% evaluate refinable integral
[x,y] = refine(mask, deriv, numRefineSteps, ...
               'FirstMaskIndex',firstMaskIndex, 'Display','on');

