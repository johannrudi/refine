%example_fct_1d_daub3  Evaluation of Daubechies' refinable function.
%  TODO
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

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

% set derivatives matrix with `num_cols = dimension`
deriv = 0;

% set number of refine steps
numRefineSteps = 5;


%% Evaluate Refinable Function

% compute refinable function
[x,y] = refine(mask, deriv, numRefineSteps);


%% Plot

% extend graph for plot
x = [-1 ; x ; 5];
y = [0 ; y ; 0];

% plot function
figure(1)
hPlot = plot(x, y);

% format plot
hTitle = title(['Daubechies'' refinable function \phi_3, ' ...
                'derivative=' num2str(deriv)]);
hXLabel = xlabel('x');
hYLabel = ylabel('y');

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

