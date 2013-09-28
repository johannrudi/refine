%EXAMPLE_FNC_1D_DAUB3 Evaluation of 1-D Daubechies' refinable function.
% Evaluates the Daubechies’ refinable function `phi_N` which generate
% orthogonal wavelets with compact support for `N = 3`.
%
% See also: EXAMPLE_FNC_1D_BSPLINE, EXAMPLE_FNC_2D, EXAMPLE_FNC_3D, REFINE.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

%% Set MATLAB Environment

% clearing
clear
clf

% initialize refine
refine_init


%% Set Parameters

% set mask coefficients
mask = [ 0.332670552950 ...
         0.806891509311 ...
         0.459877502118 ...
        -0.135011020010 ...
        -0.085441273882 ...
         0.035226291882 ...
];

% set derivative
deriv = 0;

% set number of refine steps
numRefineSteps = 5;


%% Compute Function Evaluation

% evaluate refinable function
[x,y] = refine(mask, deriv, numRefineSteps);


%% Plot

% extend data for plot
x = [-1 ; x ; 5];
y = [0 ; y ; 0];

% plot function
figure(1)
hPlot = plot(x, y);

% format plot
hTitle = title('Evaluation of 1-D Daubechies'' refinable function');
hXLabel = xlabel('x');
hYLabel = ylabel('Daubechies'' function phi_3');

set([hXLabel, hYLabel]         , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );
set( hPlot                     , ...
    'LineWidth'  , 2           );
set( gca                       , ...
    'YLim'       , [-0.5 1.5]  , ...
    'LineWidth'  , 1           );

