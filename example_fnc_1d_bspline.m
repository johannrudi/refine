%example_fct_1d_bspline  Evaluation of 1-D cardinal B-splines.
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

% set b-spline order
bsplineOrder = 4;

% set derivatives matrix with `num_cols = dimension`
deriv = 0;

% set number of refine steps
numRefineSteps = 4;

% set start index of mask coefficients
firstMaskIndex = -2;


%% Evaluate Refinable Function

% compute mask coefficients
mask = factorial(bsplineOrder) * ones(1, bsplineOrder+1) ...
       ./ factorial(bsplineOrder:-1:0) ...
       ./ factorial(0:bsplineOrder) / 2^(bsplineOrder-1);

% compute refinable function
[x,y] = refine(mask, deriv, numRefineSteps, ...
               'FirstMaskIndex',firstMaskIndex, 'Display','on');


%% Plot

% extend graph for plot
x = [x(1) - abs(x(end) - x(1))/8 ; x ; x(end) + abs(x(end) - x(1))/8];
y = [0 ; y ; 0];

% plot function
figure(1)
hPlot = plot(x(2:end-1),y(2:end-1),'or', x,y,'-b');

% format plot
hTitle = title(['B-Spline of order ' num2str(bsplineOrder) ...
                ', derivative=' num2str(deriv)]);
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
    'YLim'       , [-0.1 0.8]  , ...
    'LineWidth'  , 1           );

