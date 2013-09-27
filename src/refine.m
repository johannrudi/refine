%REFINE Evaluate refinable functions and integrals.
% Computes function values and derivatives of refinable functions as well as
% integrals of products of refinable functions with derivatives.
% The evaluation of a function or integral is based on solving an eigenvector
% problem derived from the refinement equations. The theory was developed
% in [1].
%
% A refinable function `f` satisfies the refinement equation
%
%   f(x) = \sum_{k} a_k * f(2 * x - k),
%
% with some mask coefficients `a_k`. The variables `x` and `k` are scalars or
% vectors of the same dimension `d`, where `x` is a real scalar/vector
% and `k` is an integer scalar/vector. Integrals `I(z)` of products of
% refinable functions `f_i` with derivatives `D^(m_i)` have the form
%
%   I(z) = \int f_0(x) * [D^(m_1) f_1](x - z_1)
%                      * ...
%                      * [D^(m_s) f_l](x - z_l) \dx,
%
% where `z` is an `s`-dimensional vector of integers and `s = l * d`.
% Note that `f_0` has no derivatives, so `m_0 = 0`. `I(z)` also satisfies
% a refinement equation with mask coefficients that can be computed using
% the mask coefficients of `f_i`, `i = 1,...,s`. These integrals are called
% refinable integrals.
%
% Syntax:
% [x,y] = REFINE(mask, [derivatives, numRefineSteps, ...])
%
% Input:
% mask{}[]        cell-array with mask coefficients of `f_i`, `i = 0,...,s`
% derivatives[]   derivatives `m_i` of `f_i`, `i = 1,...,s` (optional)
% numRefineSteps  number of additional refinement steps (optional)
%
% Optional Input as ('keyword', value) Pairs:
% 'FirstMaskIndex', firstMaskIndex[] (default [0,...,0])
%     sets and array with start indices of mask coefficients that determine
%     the location of the nodes at which the function/integral is evaluated
% 'EvaluateBoundary', 'on' / 'off' (default 'off')
%     forces the evaluation at nodes at the boundary (where the
%     function/integral is zero)
% 'Diplay', 'on' / 'off' (default 'off')
%     prints results of the computation
% 'WriteFile', filename
%     writes results of the computation into a file
% 'FormatArg', print_format_string
%     set format for display and file output of node values
% 'FormatVal', print_format_string (default '%+1.12e')
%     set format for display and file output of function/integral values
%
% Output:
% x[]             node values
% y[]             function/integral values
%
% Formal Definition of the Mask Coefficients:
% The mask coefficients for a refinable function need to be described by
% a matrix that is of an order that equals the dimension of the function domain.
% Furthermore, the matrix has to have the same number of elements in each
% dimension. In order to define a function
%
%   f : R^d --> R
%
% by mask coefficients stored in the matrix `M`, it is required that `M` is
% of order `d`, thus
%
%   size(M) = [1, n]                             ,   if d = 1,
%   size(M) = [n, n, ..., n], (n appears d times),   if d > 1,
%
% where `n` is the number of mask coefficients in each dimension.
%
% Formal Definition of the Derivatives and First Indices for the Mask:
% Both, the `derivatives` matrix and the `firstMaskIndex` matrix, are of the
% same size:
%
%   size(derivatives) = size(firstMaskIndex) = [1, d]  ,    if function eval.,
%   size(derivatives) = size(firstMaskIndex) = [l+1, d],    if integral eval.,
%
% where `l+1` is the number of functions in the integrand and `d` is the
% dimension of each function's domain.
%
% References:
% [1] W. Dahmen, C. A. Micchelli, "Using the refinement equation for computing
%     integrals of wavelets," SIAM J. Numer. Anal. 30, 1993, pp. 507-537.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function [x,y] = refine(mask, derivatives, numRefineSteps, varargin)


%% Check Input & Set Basic Parameters

% check number of functions constituted by mask coefficients `mask`
if isa(mask, 'cell')
    % set number of functions
    numFncs = length(mask);
    if numFncs == 1
        mask = mask{1};
    end
else
    % set number of functions to one since a matrix was given
    numFncs = 1;
end

% set dimension of the node space `dim` to the dimension of first function
if numFncs == 1 % if function evaluation
    dim = sum(size(mask) > 1);
    if isvector(mask) % if mask stored as a vector
        dim = 1;
    end
else % if integral evaluation
    dim = sum(size(mask{1}) > 1);
    if isvector(mask{1}) % if mask stored as a vector
        dim = 1;
    end
end

% check dimensions of mask coefficients of 2nd, 3rd,... function
for j = 2:numFncs % loop over all but first function
    if sum(size(mask{j}) > 1) ~= dim && ~( dim == 1 && isvector(mask{j}) )
        error(['Invalid input: Dimensions of mask ' num2str(j) ...
               ' does not match to first mask'])
    end
end

% set dimension of multi-index
if numFncs == 1 % if function evaluation
    index_dim = dim;
else % if integral evaluation
    index_dim = dim * (numFncs - 1);
end

% check number of mask coefficients
numCoeff = zeros(numFncs, 1);
if numFncs == 1 % if function evaluation
    numCoeff = length(mask);
    if dim > 1 && min(size(mask) == numCoeff) == 0
        error(['Invalid input: Number of mask coefficients ' ...
               'not equal in every dimension'])
    end
else % if integral evaluation
    for j = 1:numFncs
        numCoeff(j) = length(mask{j});
        if dim > 1 && min(size(mask{j}) == numCoeff(j)) == 0
            error(['Invalid input: Number of coefficients in mask ' ...
                   num2str(j) ' not equal in every dimension'])
        end
    end
end

% check derivatives
if nargin < 2 || ( isscalar(derivatives) && derivatives == 0 ) % if no deriv.
    derivatives = zeros(numFncs, dim);
else
    if ~isequal(size(derivatives), [numFncs dim])
        error('Invalid input: Derivatives do not match to mask coefficients')
    end
    if numFncs > 1
        % overwrite the derivatives for the first function to zero
        derivatives(1,:) = 0;
    end
end

% check number of refine steps `numRefineSteps`
if nargin < 3 || ~isscalar(numRefineSteps)
    numRefineSteps = 0;
else
    numRefineSteps = floor(numRefineSteps);
end


%% Set Preferences Passed Via `varargin`

% set default preferences
pref_firstMaskIndex = 0;
firstMaskIndex = zeros(numFncs, dim);

pref_evalBoundary = 0;

pref_display = 0;
pref_writeFile = 0;

pref_formatArg = 0;
format_arg = '%+1.12e';
pref_formatVal = 0;
format_val = '%+1.12e';

% check for keywords in input `varargin` & assign some values to parameters
for i = 1:2:length(varargin) % loop over all keywords
    % first index of mask coefficients
    if strcmpi(varargin{i}, 'FirstMaskIndex')
        if isequal(size(varargin{i+1}), [numFncs dim])
            pref_firstMaskIndex = 1;
            firstMaskIndex = varargin{i+1};
        else
            error('Invalid input: Start indices of mask coefficients')
        end

    % evaluate boundary
    elseif strcmpi(varargin{i}, 'EvaluateBoundary')
        if strcmpi(varargin{i+1}, 'on') || ...
           strcmpi(varargin{i+1}, 'true') || ...
           strcmpi(varargin{i+1}, 'yes') || ...
           isequal(varargin{i+1} == 1, 1)
            pref_evalBoundary = 1;
        end

    % display
    elseif strcmpi(varargin{i}, 'Display')
        if strcmpi(varargin{i+1}, 'on') || ...
           strcmpi(varargin{i+1}, 'true') || ...
           strcmpi(varargin{i+1}, 'yes') || ...
           isequal(varargin{i+1} == 1, 1)
            pref_display = 1;
        end

    % write to file
    elseif strcmpi(varargin{i}, 'WriteFile')
        if isa(varargin{i+1}, 'char')
            pref_writeFile = i + 1;
        else
            error('Invalid input: Filename for output file')
        end

    % format for display and file output of node values
    elseif strcmpi(varargin{i}, 'FormatArg')
        if isa(varargin{i+1}, 'char')
            pref_formatArg = i + 1;
        else
            error('Invalid input: Format of function/integral argument')
        end

    % format for display and file output of function/integral values
    elseif strcmpi(varargin{i}, 'FormatVal')
        if isa(varargin{i+1}, 'char')
            pref_formatVal = i + 1;
        else
            error('Invalid input: Format of function/integral value')
        end

    % unknown keyword
    else
        error(['Invalid input: Unknown keyword ' varargin{i}])
    end % of `varargin` parsing
end

% assign values of input `varargin` to parameters
if pref_formatArg
    % set format of node values
    format_arg = varargin{pref_formatArg};
end
if pref_formatVal
    % set format of function/integral values
    format_val = varargin{pref_formatVal};
end
if pref_display
    % start measuring time
    tic

    % set shorter format of the argument
    if ~pref_formatArg % if format of argument not customized
        if numRefineSteps == 0 % if evaluation at integer nodes only
            format_arg = '%+2i';
        else
            format_arg = '%+2.4f';
        end
    end

    % set format of display output
    format_display = [ ...
        '(' strtrim(repmat([format_arg ' '], 1, index_dim)) ') = ' ...
        format_val ...
    ];
end
if pref_writeFile
      % set filename
      filename = varargin{pref_writeFile};

      % set format of file output
      format_writeFile = [ ...
          repmat([format_arg '  '], 1, index_dim) '=  ' format_val ...
      ];
end


%% Set Parameters for Output

if pref_display || pref_writeFile
    if numFncs == 1
        % set name of function
        disp_name = refine_print_getStringOfFunction ('f', derivatives);
    else
        % set name of integral
        disp_name = 'I';

        % set argument of the integral function
        disp_arg = repmat(' z_#', 1, numFncs-1);
        disp_arg(4:4:end) = num2str(1:numFncs-1, '%i');
    end

    % set space of nodes
    if numRefineSteps == 0
        disp_nodeSpace = ['Z^' num2str(dim)];
    else
        disp_nodeSpace = ...
            ['(Z^' num2str(dim) ')/(2^' num2str(numRefineSteps) ')'];
    end
end


%% Output On Screen

if pref_display
    if numFncs == 1 % if function evaluation
        % print input on screen
        refine_printl('INPUT: Mask coefficients of one refinable function f');
        refine_printl(['- dimension of function domain = ' num2str(dim)]);
        refine_printl(['- number of mask coefficients in each dimension = ' ...
                       num2str(numCoeff)]);
        if pref_firstMaskIndex
            refine_printl(['- start index of mask coefficients = ' ...
                           num2str(firstMaskIndex)]);
        end
        refine_printl(['- derivative(s) = ' num2str(derivatives)]);
        refine_printl(['- number of refine steps = ' num2str(numRefineSteps)]);
        refine_printl('');

        % print output on screen
        refine_printl('OUTPUT: Evaluation of function');
        refine_printl(['  ' disp_name '(x)']);
        refine_printl(['with x in ' disp_nodeSpace]);
        refine_printl('');
    else % if integral evaluation
        % print input on screen
        refine_printl(['INPUT: Mask coefficients of ' num2str(numFncs) ' ' ...
                       'refinable functions f_i, i=0,...,' num2str(numFncs-1)]);
        refine_printl(['- dimension of function domain = ' num2str(dim)]);
        refine_printl(['- f_0: number of mask coefficients in ' ...
                       'each dimension = ' num2str(numCoeff(1))]);
        for k = 2:numFncs % loop over all but first function
            refine_printl(['- f_' num2str(k-1) ...
                           ': number of mask coefficients in ' ...
                           'each dimension = ' num2str(numCoeff(k))]);
            if pref_firstMaskIndex
                refine_printl(['- f_' num2str(k-1) ...
                               ': start index of mask coefficients = ' ...
                               num2str(firstMaskIndex(k,:))]);
            end
            refine_printl(['- f_' num2str(k-1) ...
                           ': derivative(s) = ' num2str(derivatives(k,:))]);
        end
        refine_printl(['- number of refine steps = ' num2str(numRefineSteps)]);
        refine_printl('');

        % print output on screen
        refine_printl('OUTPUT: Evaluation of integral');
        disp_integrand = 'f_0(x)';
        for k = 2:numFncs % loop over all but first function
            str_fnc = refine_print_getStringOfFunction ( ...
                          ['f_' num2str(k - 1)], derivatives(k,:));
            disp_integrand = [disp_integrand ...
                              ' * ' str_fnc '(x-z_' num2str(k - 1) ')'];
        end
        refine_printl(['  ' disp_name '(' disp_arg(2:end) ') = ' ...
                       '\\int ' disp_integrand ' \\dx']);
        refine_printl(['with z_i in ' disp_nodeSpace ...
                       ', i=1,...,' num2str(numFncs - 1)]);
        refine_printl('');
    end
end


%% Set Parameters for Computation

% set multi-indices, modify derivatives and mask coefficients
if numFncs == 1 % if evaluation of function
    % normalize mask coefficients
    mask = refine_scaleMask(mask, dim, 1);

    % transform derivatives input into a column vector
    derivatives = derivatives(:);

    % set first index
    index_first = firstMaskIndex(:);

    % set last index
    index_last = firstMaskIndex(:) + numCoeff - 1;
else % if evaluation of integral
    % normalize mask coefficients
    for k = 1:numFncs
        mask{k} = refine_scaleMask(mask{k}, dim, 1);
    end

    % transform derivatives matrix to column vector
    derivatives = reshape(derivatives(2:end,:)', (numFncs - 1) * dim, 1);

    % set multi-indices for mask coefficients of all functions
    IndexMask_first = firstMaskIndex;
    IndexMask_last = firstMaskIndex + repmat(numCoeff, 1, dim) - 1;

    % set first index for mask coefficients of the integral
    index_first = repmat(IndexMask_first(1,:), 1, numFncs - 1)' - ...
                  reshape(IndexMask_last(2:end,:)', (numFncs - 1) * dim, 1);

    % set last index for the mask coefficients of the integral
    index_last = repmat(IndexMask_last(1,:), 1, numFncs - 1)' - ...
                 reshape(IndexMask_first(2:end,:)', (numFncs - 1) * dim, 1);

    % compute mask coefficients of the integral
    mask = refine_computeIntMask(dim, numFncs, mask, derivatives, ...
                                 index_first, index_last, ...
                                 IndexMask_first, IndexMask_last);
end

% set parameters depending on number of refine steps `numRefineSteps`
if numRefineSteps == 0 && ~pref_evalBoundary
    % set number of nodes (= number of x- and y-values)
    numNodes = prod(abs(index_last - index_first) - 1);

    % create multi-index for looping over all x-values
    muid = multiindex_create(index_first + 1, index_last - 1);
else
    % set number of nodes (= number of x- and y-values)
    numNodes = prod(2^numRefineSteps * abs(index_last - index_first) + 1);

    % create multi-index for looping over all x-values
    muid = multiindex_create(index_first, index_last);
end

% set value that the x-values are incremented by
incrementBy = 1 / (2^numRefineSteps);


%% Compute y-Values

if numRefineSteps == 0 && ~pref_evalBoundary
    % compute just inner y-values; faster, because y-values on boundary
    % that equal zero anyways are not computed
    y = refine_computeInitialValues(dim, numFncs, ...
                                    mask, derivatives, ...
                                    index_first, index_last);
else
    % compute inner y-values and y-values on boundary
    y = refine_computeRefinedValues(dim, numFncs, ...
                                    mask, derivatives, ...
                                    index_first, index_last, ...
                                    numRefineSteps);
end


%% Compute x-Values & Output of Computed Values

% if applicable, open file for output & write header
if pref_writeFile
    file = fopen(filename, 'w');
    fprintf(file, '/**\n');
    fprintf(file, [' * argument(s) of ' disp_name ...
                   ' = value of ' disp_name '\n']);
    %TODO maybe write mask, deriv, etc.
    fprintf(file, ' */\n\n');
end

% create empty matrix for x-values
x = zeros(numNodes, max([1, numFncs-1]) * dim);

% set all x-values, print x- and y-values
row = 1;
isEnd = 0;
while ~isEnd % loop over all x-values
    % get multi-index
    index = multiindex_getPos(muid);

    % set x-value
    x(row,:) = index;

    % if applicable, output on screen
    if pref_display && abs(y(row)) > 1000 * eps
        refine_printl([disp_name num2str([index' y(row)], format_display)]);
    end

    % if applicable, write to file
    if pref_writeFile
        fprintf(file, [format_writeFile '\n'], index, y(row));
    end

    % increment multi-index
    [muid,isEnd] = multiindex_increment_rowMajor(muid, incrementBy);

    % increment row
    row = row + 1;
end


%% Finalize

% if applicable, close file for output
if pref_writeFile
    fclose(file);
end


% stop Measuring Time
if pref_display
    refine_printl(['FINISHED. Time: ' num2str(toc) ' sec.']);
end

% end function
end
