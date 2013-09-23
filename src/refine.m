%refine TODO
%  Computes integrals of products of 1-dimensional refinable functions on
%  integer nodes.
%
%  Syntax:
%  [x,y] = refine(mask, derivatives, numRefineSteps, ...TODO)
%
%  Input:
%  mask{}[]          cell-array with mask coefficients of `phi_i`, `i=0,...,s`
%  firstMaskIndex[]  array with start indices of coefficients
%  derivatives[]     derivatives `mu_i` of `phi_i`, `i=1,...,s` (optional)
%  output            output preferences string (optional); 's': output on screen
%
%  Output:
%  X[]                x-values
%  y[]                y-values
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-06-29
%  ----------------------------------------------------------------------------

function [x,y] = refine(mask, derivatives, numRefineSteps, varargin)

%TODO rewrite display output
DISPLAY_PREFIX = '[refine]';


%% Check Input

% check number of functions constituted by mask coefficients `mask`
if isa(mask, 'cell')
    % set number of functions
    numFncs = length(mask);
    if numFncs == 1
        mask = mask{1};
    end
else
    % set number of functions
    numFncs = 1;
end

% set dimension `dim` to the dimension of first function
if numFncs == 1 % if just one function
    dim = sum(size(mask) > 1);
    if isvector(mask) % if mask stored in vector
        dim = 1;
    end
else
    dim = sum(size(mask{1}) > 1);
    if isvector(mask{1}) % if mask stored in vector
        dim = 1;
    end
end

% check dimensions of mask coefficients of 2nd, 3rd,... function
for j = 2:numFncs % loop over all but first function
    if sum(size(mask{j}) > 1) ~= dim && ~(dim == 1 && isvector(mask{j}))
        error(['Invalid input: Dimensions of mask ' ...
               num2str(j) ' does not match.'])
    end
end

% set dimension of multi-index
if numFncs == 1 % if just one function
    index_dim = dim;
else
    index_dim = dim * (numFncs - 1);
end

% check number of mask coefficients
numCoeff = zeros(numFncs, 1);
if numFncs == 1
    numCoeff = length(mask);
    if dim > 1 && min(size(mask) == numCoeff) == 0
        error(['Invalid input: Number of mask coefficients ' ...
               'not equal in every dimension.'])
    end
else
    for j = 1:numFncs
        numCoeff(j) = length(mask{j});
        if dim > 1 && min(size(mask{j}) == numCoeff(j)) == 0
            error(['Invalid input: Number of coefficients in mask ' ...
                   num2str(j) ' not equal in every dimension.'])
        end
    end
end

% check derivatives `derivatives`
if nargin < 2 || (isscalar(derivatives) && derivatives == 0) % if no deriv.
    derivatives = zeros(numFncs, dim);
else
    if ~isequal(size(derivatives), [numFncs dim])
        error('Invalid input: Derivatives do not match to mask coefficients.')
    end
    if numFncs > 1
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
for i = 1:2:length(varargin) % loop over all property keywords
    % first index of mask coefficients
    if strcmpi(varargin{i}, 'FirstMaskIndex')
        if isequal(size(varargin{i+1}), [numFncs dim])
            pref_firstMaskIndex = 1;
            firstMaskIndex = varargin{i+1};
        else
            error('Invalid input: Start indices of mask coefficients.')
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
            error('Invalid input: Filename for output file.')
        end

    % format of function/integral argument used for display and file output
    elseif strcmpi(varargin{i}, 'FormatArg')
        if isa(varargin{i+1}, 'char')
            pref_formatArg = i + 1;
        else
            error('Invalid input: Format of function/integral argument.')
        end

    % format of function/integral value used for display and file output
    elseif strcmpi(varargin{i}, 'FormatVal')
        if isa(varargin{i+1}, 'char')
            pref_formatVal = i + 1;
        else
            error('Invalid input: Format of function/integral value.')
        end

    %TODO else unknown argument
    end % of `varargin` check
end

% assign remaining values of input `varargin` to parameters
if pref_formatArg
    % set format of function/integral argument
    format_arg = varargin{pref_formatArg};
end
if pref_formatVal
    % set format of function/integral value
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
        % write derivatives to string
        str_deriv = num2str(derivatives, ',%d');

        % set name of function
        if isequal(derivatives, zeros(numFncs, dim)) % if no derivatives
            disp_name = 'phi';
        else
            if dim == 1 % if 1-dimensional
                disp_name = ['(D^' str_deriv(2:end) ' phi)'];
            else
                disp_name = ['(D^(' str_deriv(2:end) ') phi)'];
            end
        end
    else
        % set name of integral function
        disp_name = 'H';
    end
end


%% Output On Screen

if pref_display
    if numFncs == 1
        % write derivatives to string
        str_deriv = num2str(derivatives, ',%d');

        % set name of function
        if isequal(derivatives, zeros(numFncs, dim)) % if no derivatives
            disp_name = 'phi';
        else
            if dim == 1 % if 1-dimensional
                disp_name = ['(D^' str_deriv(2:end) ' phi)'];
            else
                disp_name = ['(D^(' str_deriv(2:end) ') phi)'];
            end
        end

        % print input on screen
        disp([DISPLAY_PREFIX ' INPUT: Mask coefficients of one ' ...
              'refinable function phi'])
        disp([DISPLAY_PREFIX '   - dimension = ' num2str(dim)])
        disp([DISPLAY_PREFIX '   - number of mask coefficients in each ' ...
              'dimension = ' num2str(numCoeff)])
        if pref_firstMaskIndex
            disp([DISPLAY_PREFIX '   - start index of mask coefficients = ' ...
                  num2str(firstMaskIndex)])
        end
        disp([DISPLAY_PREFIX '   - derivative(s) = ' num2str(derivatives)])
        disp([DISPLAY_PREFIX '   - number of refine steps = ' ...
              num2str(numRefineSteps)])
        disp(DISPLAY_PREFIX)

        % print output on screen
        disp([DISPLAY_PREFIX ' OUTPUT: Evaluation of function'])

        % print function on screen
        disp([DISPLAY_PREFIX '     ' disp_name '(x)'])
        fprintf([DISPLAY_PREFIX ' with x element of '])
        if numRefineSteps == 0
            fprintf('Z^%i', dim)
        else
            fprintf('(Z^%i)/(2^%i)', dim, numRefineSteps)
        end
        fprintf('\n[refine]\n')
    else
        % write derivatives to string
        str_deriv = num2str(derivatives, ',%1d');

        % write argument of integral function
        str_arg = repmat(' z_#', 1, numFncs-1);
        str_arg(4:4:end) = num2str(1:numFncs-1, '%i');

        % set name of integral function
        disp_name = 'H';

        % print input on screen
        disp([DISPLAY_PREFIX ' INPUT: Mask coefficients of ' num2str(numFncs) ...
              ' refinable functions phi_i, i=0,...,' num2str(numFncs-1)])
        disp([DISPLAY_PREFIX '   - dimension = ' num2str(dim)])
        disp([DISPLAY_PREFIX '   - phi_0: ' ...
              'number of mask coefficients in each dimension = ' ...
              num2str(numCoeff(1))])

        for k = 2:numFncs % loop over all but first function
            disp([DISPLAY_PREFIX '   - phi_' num2str(k-1) ...
                  ': number of mask coefficients in each dimension = ' ...
                  num2str(numCoeff(k))])
            if pref_firstMaskIndex
                disp([DISPLAY_PREFIX '   - phi_' num2str(k-1) ...
                      ': start index of mask coefficients = ' ...
                      num2str(firstMaskIndex(k,:))])
            end
            disp([DISPLAY_PREFIX '   - phi_' num2str(k-1) ...
                  ': derivative(s) = ' num2str(derivatives(k,:))])
        end
        disp([DISPLAY_PREFIX '   - number of refine steps = ' ...
              num2str(numRefineSteps)])
        disp(DISPLAY_PREFIX)

        % print output on screen
        disp([DISPLAY_PREFIX ' OUTPUT: Evaluation of integral'])

        % print integral function on screen
        fprintf([DISPLAY_PREFIX '     ' disp_name '(' str_arg(2:end) ...
                 ') = INT phi_0(x)'])
        for k = 2:numFncs % loop over all but first function
            if dim == 1
                fprintf([' * (D^' str_deriv(k,2:end) ' phi_%i)(x-z_%i)'], ...
                        k-1, k-1)
            else
                fprintf([' * (D^(' str_deriv(k,2:end) ') phi_%i)(x-z_%i)'],...
                        k-1, k-1)
            end
        end
        fprintf(' dx\n')
        fprintf([DISPLAY_PREFIX ' with x element of R^%i, ', dim])
        if numRefineSteps == 0
            fprintf('and z_i element of Z^%i, ', dim)
        else
            fprintf('and z_i element of (Z^%i)/(2^%i), ', dim, numRefineSteps)
        end
        fprintf('for all i=1,...,%i\n[refine]\n', numFncs-1)
    end
end


%% Set Parameters for Computation

% set multi-indices, modify derivatives and mask coefficients
if numFncs == 1 % if function values should be computed
    % norm mask coefficients
    mask = refine_scaleMask(mask, dim, 1);

    % transform derivatives input into a column vector
    derivatives = derivatives(:);

    % set first index
    index_first = firstMaskIndex(:);

    % set last index
    index_last = firstMaskIndex(:) + numCoeff - 1;
else % if evaluation of integral should be computed
    % norm mask coefficients
    for k = 1:numFncs
        mask{k} = refine_scaleMask(mask{k}, dim, 1);
    end

    % transform derivatives matrix `derivatives` to column vector
    derivatives = reshape(derivatives(2:end,:)', (numFncs - 1) * dim, 1);

    % set multi-indices for mask coefficients of all functions
    IndexMask_first = firstMaskIndex;
    IndexMask_last = firstMaskIndex + repmat(numCoeff, 1, dim) - 1;

    % set first index for mask coefficients of the integral TODO
    index_first = repmat(IndexMask_first(1,:), 1, numFncs - 1)' - ...
                  reshape(IndexMask_last(2:end,:)', (numFncs - 1) * dim, 1);

    % set last index for the mask coefficients of the integral TODO
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
    % compute just inner y-values; faster, because y-values on boundary,
    % that equal zero anyways, are not computed
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


%% Compute x-Values

% if applicable, open file for output & write header
if pref_writeFile
    file = fopen(filename, 'w');
    fprintf(file, '/**\n');
    fprintf(file, [' * argument(s) of ' disp_name ...
                   ' = value of ' disp_name '\n']);
    %TODO write mask, deriv, etc.
    fprintf(file, ' */\n\n');
end

% create empty matrix for x-values
x = zeros(numNodes, max([1, numFncs-1]) * dim);

% set all x-values & print x- and y-values
row = 1;
isEnd = 0;
while ~isEnd % loop over all x-values
    % get multi-index
    index = multiindex_getPos(muid);

    % set x-value
    x(row,:) = index;

    % if applicable, output on screen
    if pref_display && abs(y(row)) > 1000 * eps
        fprintf([DISPLAY_PREFIX ' ' disp_name format_display '\n'], ...
                index, y(row));
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

% if applicable, close file for output
if pref_writeFile
    fclose(file);
end


%% Stop Measuring Time
if pref_display
    disp([DISPLAY_PREFIX ' FINISHED. Time: ' num2str(toc) ' sec.'])
end

% end function
end
