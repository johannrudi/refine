%refine_computeInitialValues  Compute initial values.
%  TODO
%  Computes integrals of products of 1-dimensional refinable functions on
%  integer nodes.
%
%  Syntax:
%  [X,y] = refine_computeInitialValues(mask, maskIndexStart, derivatives, output)
%
%  Input:
%  mask{}[]          cell-array with mask coefficients of `phi_i`, `i=0,...,s`
%  maskIndexStart[]  array with start indices of coefficients
%  derivatives[]     derivatives `mu_i` of `phi_i`, `i=1,...,s` (optional)
%  output            output preferences string (optional); 's': output on screen
%
%  Output:
%  X[]                x-values
%  y[]                y-values
%
%  Description:
%  A refinable function `phi` satisfies a refinement equation
%    phi(x) = \sum_{k \in ZZ} a_j phi(2 * x - k),   x \in RR,
%  with mask coefficients `a_j`.
%
%  Integrals `H(alpha)` of products of refinable functions `phi_i`
%    H(alpha) = \int phi_0(x) * (D^{mu_1} phi_1)(x - alpha_1) * ...
%                   * (D^{mu_s} phi_s)(x - alpha_s) dx,   alpha \in ZZ^s,
%  also satisfy a refinement equation with mask coefficients that can be
%  computet using the mask coefficients of `phi_i`, `i=1,...,s`. These
%  integrals are called refinable integrals.
%
%  References:
%  [1] W. Dahmen, C. A. Micchelli, "Using the refinement equation for computing
%      integrals of wavelets," SIAM J. Numer. Anal. 30, 1993, pp. 507-537.
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-06-29
%  ----------------------------------------------------------------------------

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
% `absNumDeriv` using the formula: TODO: does not seem to be what's computed
%   #possibilities to choose `index_dim` out of `n` non-distinguishable objects
n = 0:absNumDeriv;
numPoss = sum( ...
              factorial(n + (index_dim - 1) * ones(size(n))) ...
              ./ factorial(n) ...
              ./ factorial((index_dim - 1) * ones(size(n))) ...
);


%% Create System of Linear Equations with Matrix `A` and Right-Hand Side `b`

% create empty matrix `A` and right-hand side `b`
A = zeros(index_numElements + numPoss, index_numElements);
b = zeros(index_numElements + numPoss, 1);

% create multi-index for looping over matrix rows or columns
%TODO delete this? muid_row = multiindex_create(index_first, index_last);
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

%TODO maybe use diag?
% subtract derivatives constant `derivConst` from diagonal
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
% gives us y-values on integer x-values

% solve
y = A \ b;

if numFncs > 1 % if values of integral are computed
    % scale solution of `y` at integer nodes according to 3rd equation in [1]
    % on p. 528
    y = y * (-1)^absNumDeriv;
end

% end function
end
