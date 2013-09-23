%example_fct_2d  Evaluation of 2-D box spline.
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
mask = [
    0.125  0.25   0.125   0.0   ; ...
    0.25   0.625  0.5     0.125 ; ...
    0.125  0.5    0.625   0.25  ; ...
    0.0    0.125  0.25    0.125 ...
];

% set derivatives vector (number of columns = dimension of function)
deriv = [1 0];

% set number of refine steps
numRefineSteps = 3;

% set first index of mask coefficients
% (number of columns = dimension of function)
firstMaskIndex = [-1 -2];


%% Evaluate Refinable Function

% compute refinable function
[x,y] = refine(mask, deriv, numRefineSteps, 'FirstMaskIndex',firstMaskIndex);


%% Plot

% prepare data
numKnots_oneDim = (length(mask) - 1) * 2^numRefineSteps + 1;
X = reshape(x(:,1), numKnots_oneDim, numKnots_oneDim);
Y = flipud(reshape(x(:,2), numKnots_oneDim, numKnots_oneDim));
Z = reshape(y, numKnots_oneDim, numKnots_oneDim);

% plot function
figure(1)
hPlot = surf(X, Y, Z);

% format plot
hTitle = title(['Bivariate box spline, derivatives=' num2str(deriv, ' %i')]);
hXLabel = xlabel('x');
hYLabel = ylabel('y');
hZLabel = zlabel('z');

set([hXLabel, hYLabel]         , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

