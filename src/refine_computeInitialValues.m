%REFINE_COMPUTEINITIALVALUES Evaluate refinable fnc/integral at integer nodes.
% Computes function values and derivatives of refinable functions at
% integer nodes as well as integrals of products of refinable functions with
% derivatives at integer nodes. The evaluation of a function or integral is
% based on solving an eigenvector problem derived from the refinement equations.
%
% Syntax:
% y = REFINE_COMPUTEINITIALVALUES(dim, numFncs, mask, derivatives, ...
%                                 index_first, index_last)
%
% Input:
% dim            dimension of the functions' domains
% numFncs        number of functions
% mask{}[]       cell-array with mask coefficients of the functions
% derivatives[]  derivatives of the functions
% index_first[]  beginning of multi-index
% index_last[]   end of multi-index
%
% Output:
% y[]            function/integral values
%
% See also: REFINE_COMPUTEREFINEDVALUES, REFINE.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function y = refine_computeInitialValues(dim, numFncs, ...
                                         mask, derivatives, ...
                                         index_first, index_last)

%% Set Parameters

% set index dimension
index_dim = length(index_first);

% set number of index elements
index_numElements = prod(abs(index_last - index_first - 2) + 1);

% set total number of all derivatives
absNumDeriv = sum(sum(abs(derivatives)));

% set constant depending on derivatives
derivConst = 1 / (2^absNumDeriv);

% compute the number of possibilities of derivatives lower or equal
n = 0:absNumDeriv;
numPoss = sum( ...
    factorial(n + (index_dim - 1) * ones(size(n))) ./ factorial(n) ...
    ./ factorial((index_dim - 1) * ones(size(n))) ...
);


%% Create System of Linear Equations with Matrix `A` and Right-Hand Side `b`

% create empty matrix `A` and right-hand side `b`
A = zeros(index_numElements + numPoss, index_numElements);
b = zeros(index_numElements + numPoss, 1);

% create multi-index for looping over matrix rows or columns
muid_row = multiindex_create(index_first + 1, index_last - 1);

% set upper part of matrix `A`
row = 1;
isEnd_row = 0;
while ~isEnd_row % loop over all row multi-indices
    % get row multi-index
    index_row = multiindex_getPos(muid_row);

    % create multi-indices for looping over matrix columns
    muid_col = multiindex_create( ...
                   max(index_first + 1, 2 * index_row - index_last), ...
                   min(index_last - 1, 2 * index_row - index_first));

    isEnd_col = 0;
    while ~isEnd_col % loop over all column multi-indices
        % get column multi-index
        index_col = multiindex_getPos(muid_col);

        % get matrix column
        col = multiindex_nDimTo1Dim_rowMajor(index_col, ...
                                             index_first + 1, ...
                                             index_last - 1) + 1;

        % get 1-dimensional mask index
        index = 2 * index_row - index_col;
        if numFncs == 1
            index_1dim = multiindex_nDimTo1Dim_colMajor(index, ...
                                                        index_first, ...
                                                        index_last) + 1;
        else
            index_1dim = multiindex_nDimTo1Dim_rowMajor(index, ...
                                                        index_first, ...
                                                        index_last) + 1;
        end

        % set matrix entry
        A(row,col) = mask(index_1dim);

        % increment column multi-index
        [muid_col,isEnd_col] = multiindex_increment_rowMajor(muid_col);
    end

    % increment row multi-index
    [muid_row,isEnd_row] = multiindex_increment_rowMajor(muid_row);

    % increment matrix row
    row = row + 1;
end

% subtract derivatives constant `derivConst` from diagonal
%TODO maybe use diag()?
diag_selector = 1 : index_numElements + numPoss + 1 : ...
                index_numElements * (index_numElements + numPoss);
A(diag_selector) = A(diag_selector) - derivConst;

% create multi-index for looping over derivatives
muid_deriv = boundedMultiindex_create(absNumDeriv, ...
                                      zeros(max([1, numFncs-1]) * dim, 1));

% set lower part of matrix `A` and right-hand side `b`
row = index_numElements + 1;
while ~boundedMultiindex_isEnd(muid_deriv) % loop over all derivatives
                                           % lower or equal than `absNumDeriv`
    % get derivatives multi-index, which determines row
    index_deriv = boundedMultiindex_getPos(muid_deriv);

    % set column multi-index (using previous row multi-index)
    muid_col = multiindex_setPosToFirst(muid_row);

    col = 1;
    isEnd_col = 0;
    while ~isEnd_col % loop over all column multi-indices
        % get column multi-index
        index_col = multiindex_getPos(muid_col);

        % set matrix entry
        A(row,col) = prod((-index_col) .^ index_deriv);

        % increment column multi-index
        [muid_col,isEnd_col] = multiindex_increment_rowMajor(muid_col);

        % increment matrix column
        col = col + 1;
    end

    % set right-hand side entry
    if isequal(index_deriv, derivatives) % if derivative multi-index and
                                         % derivative vector are equal
        b(row) = prod(factorial(derivatives));
    end

    % increment derivative multi-index
    muid_deriv = boundedMultiindex_increment(muid_deriv);

    % increment matrix row
    row = row + 1;
end


%% Solve System of Linear Equations

% solve; gives us y-values on integer x-values
y = A \ b;

if numFncs > 1 % if integral evaluation
    % scale solution of `y` at integer nodes according to 3rd equation in
    % [Dahmen, Micchelli, 1993, p. 528]
    y = y * (-1)^absNumDeriv;
end

% end function
end
