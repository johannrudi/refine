%EXAMPLE_FNC_1D_BSPLINE Evaluation of 1-D cardinal B-Spline function.
% Evaluates a univariate cardinal B-Spline, `N_k(x)`, of order `k`.
% These are B-Splines on a grid with uniform spacing, they satisfy the
% refinement relation
%
%   N_k(x) = 2^{1-k} * \sum_{j=0..k} (k choose j) N_k(2*x - j),
%
% where x is a real scalar. Hence, the mask is the (row) vector
%
%   M = 2^{k-1} * [(k choose 0), (k choose 1), ..., (k choose k)].
%
% See also: EXAMPLE_FNC_1D_DAUB3, EXAMPLE_FNC_2D, EXAMPLE_FNC_3D, REFINE.
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

% set B-Spline order
bsplineOrder = 4;

% set derivative of B-Spline
deriv = 0;

% set number of recursive refine steps
numRefineSteps = 4;

% set start index for mask coefficients s.t. B-Spline is centered about zero
firstMaskIndex = -2;


%% Compute Function Evaluation

% compute mask coefficients
mask = factorial(bsplineOrder) * ones(1, bsplineOrder + 1) ...
       ./ factorial(bsplineOrder:-1:0) ...
       ./ factorial(0:bsplineOrder) / 2^(bsplineOrder - 1);

% evaluate refinable function
[x,y] = refine(mask, deriv, numRefineSteps, ...
               'FirstMaskIndex',firstMaskIndex, 'Display','on');


%% Plot

% extend data for plot
x = [x(1) - abs(x(end) - x(1))/8 ; x ; x(end) + abs(x(end) - x(1))/8];
y = [y(1) ; y ; y(end)];

% plot function
figure(1)
hPlot = plot(x(2:end-1),y(2:end-1),'or', x,y,'-b');

% format plot
hTitle = title('Evaluation of refinable 1-D cardinal B-Spline function');
hXLabel = xlabel('x');
hYLabel = ylabel(['B-Spline, order ' num2str(bsplineOrder) ...
                  ', derivative ' num2str(deriv)]);

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

