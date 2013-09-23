%refine_computeRefinedValues  Compute refined values recursively.
%  TODO
%
%  Syntax:
%  y = refine_computeRefinedValues(dim, numFncs,
%                                  mask, derivatives,
%                                  index_first, index_last,
%                                  numRefineSteps)
%
%  Input:
%  mask{}[]          cell-array with mask coefficients of `phi_i`, `i=0,...,s`
%  maskIndexStart[]  array with start indices of coefficients
%  derivatives[]     derivatives `mu_i` of `phi_i`, `i=1,...,s` (optional)
%  output            output preferences string (optional); 's': output on screen
%
%  Output:
%  y[]                y-values
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

function y = refine_computeRefinedValues(dim, numFncs, ...
                                         mask, derivatives, ...
                                         index_first, index_last, ...
                                         numRefineSteps)

%% Set Parameters

% set number of nodes
numNodes = prod(abs(index_last - index_first) + 1);


%% If `numRefineSteps == 0`: Compute Initial y-Values On Integer x-Values

if numRefineSteps == 0
    % create vector for y-values (including y-values on border)
    y = zeros(numNodes, 1);

    % compute initial (inner, i.e. not on border) y-values on integer nodes
    y_init = refine_computeInitialValues(dim, numFncs, ...
                                         mask, derivatives, ...
                                         index_first, index_last);

    % create multi-index for looping over `y_init`
    muid = multiindex_create(index_first + 1, index_last - 1);

    % write y-values from `y_init` into `y`
    row_init = 1;
    isEnd = 0;
    while ~isEnd % loop over all multi-indices pertaining to the initial values
        % get multi-index
        index_init = multiindex_getPos(muid);

        % get row for `y`
        row = multiindex_nDimTo1Dim_rowMajor(index_init, ...
                                             index_first, ...
                                             index_last) + 1;

        % write y-value
        y(row) = y_init(row_init);

        % increment multi-index
        [muid,isEnd] = multiindex_increment_rowMajor(muid);

        % increment `row_init`
        row_init = row_init + 1;
    end

    return
end


%% Set Additional Parameters

% create index of zeros
index_zeros = zeros(length(index_first), 1);

% set total number of all derivatives & define a factor used in computation
% below
absNumDeriv = sum(sum(abs(derivatives)));
derivFactor = 2^absNumDeriv;

% set number of nodes (= number of x- and y-values)
numNodes = prod(2^numRefineSteps * abs(index_last - index_first) + 1);

% set value that the y-values are incremented by
incrementBy = 1 / (2^numRefineSteps);


%% If `numRefineSteps > 0`: Refine Recursively

% create vector for y-values
y = zeros(numNodes, 1);

% compute y-values on coarse grid with `numRefineSteps - 1` refine steps
y_coarse = refine_computeRefinedValues(dim, numFncs, ...
                                       mask, derivatives, ...
                                       index_first, index_last, ...
                                       numRefineSteps - 1);

% create multi-index for looping over `y`
muid = multiindex_create(index_first, index_last);

% write y-values from `y_coarse` into `y` & compute y-values on fine grid
row = 1;
isEnd = 0;
while ~isEnd % loop over all multi-indices
    % get current multi-index
    index = multiindex_getPos(muid);

    if isequal(mod(index, 2 * incrementBy), index_zeros) % if index on
                                                         % coarse grid
        % get row for `y_coarse`
        row_coarse = multiindex_nDimTo1Dim_incrByPowOfTwo( ...
                     	index, ...
                        index_first, ...
                        index_last, ...
                        2 * incrementBy) + 1;

        % write coarse y-value
        y(row) = y_coarse(row_coarse);
    else
        % set refinement multi-index for computing y-value on fine grid with
        % refinement equation
        index_ref_first = max(index_first, ceil(2 * index - index_last));
        index_ref_last = min(index_last, floor(2 * index - index_first));
        muid_ref = multiindex_create(index_ref_first, index_ref_last);

        isEnd_ref = 0;
        while ~isEnd_ref % loop over all refinement multi-indices
            % get refinement multi-index
            index_ref = multiindex_getPos(muid_ref);

            % get row for `y_coarse`
            row_coarse = multiindex_nDimTo1Dim_incrByPowOfTwo( ...
                               2 * index - index_ref, ...
                               index_first, ...
                               index_last, ...
                               2 * incrementBy) + 1;

            % get 1-dimensional index for mask coefficient
            if numFncs == 1
                index_mask_1dim = multiindex_nDimTo1Dim_colMajor( ...
                                      index_ref, ...
                                      index_first, ...
                                      index_last) + 1;
            else
                index_mask_1dim = multiindex_nDimTo1Dim_rowMajor( ...
                                      index_ref, ...
                                      index_first, ...
                                      index_last) + 1;
            end

            % add to `y`
            y(row) = y(row) + mask(index_mask_1dim) * y_coarse(row_coarse);

            % increment refinement multi-index
            [muid_ref,isEnd_ref] = multiindex_increment_rowMajor(muid_ref);
        end

        % scale by factor depending on derivatives
        y(row) = y(row) * derivFactor;
    end

    % increment multi-index
    [muid,isEnd] = multiindex_increment_rowMajor(muid, incrementBy);

    % increment `row`
    row = row + 1;
end

% end function
end
