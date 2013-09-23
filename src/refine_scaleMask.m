%refine_scaleMask  Scale mask coefficients.
%  Scales the mask coefficients of a refinable function to a desired value.
%  The scaling is done according to the theory in [1] such that
%    sum_{b} c_{a - 2*b}
%  attains the desired value for all integers `a` (in the univariate case),
%  where `c_{a - 2*b}` is a mask coefficient and `b` is an integer.
%  In the multivariate case integers are replaced by multi-dimensional integers.
%
%  [1] W. Dahmen, C. A. Micchelli, "Using the refinement equation for computing
%      integrals of wavelets," SIAM J. Numer. Anal. 30, 1993, pp. 507-537.
%
%  Syntax:
%  mask = refine_scaleMask(mask, dim, scaleTo)
%
%  Input:
%  mask[]   mask coefficients of a refinable function
%  dim      dimension of domain of refinable function
%  scaleTo  value that the mask coefficients will be scaled to
%
%  Output:
%  mask[]   scaled mask coefficients
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

function mask = refine_scaleMask(mask, dim, scaleTo)

%% Scaling of Univariate Functions

if dim == 1
    for k = 1:2
        % compute sum `s`
        s = sum(mask(k:2:end));

        % scale if necessary
        if abs(s - scaleTo) > 10 * eps
            mask(k:2:end) = scaleTo * mask(k:2:end) / s;
        end
    end

    return
end


%% Scaling of Multivariate Functions

% number of mask coefficients in one dimension
numCoeff = length(mask);

% set up multi-index for offset
index_offset_first = ones(dim, 1);
index_offset_last = 2 * ones(dim, 1);

% create multi-index for offset, which will be the first index when looping
% over mask coefficients
muid_offset = multiindex_create(index_offset_first, index_offset_last);

% set variables for multi-index, which will be used for looping over mask
% coefficients
index_first = ones(dim, 1);
index_last = numCoeff * ones(dim, 1);
index_length = abs(index_last - index_first) + 1;

isEnd_offset = 0;
while ~isEnd_offset % loop over all multi-indices for offset
    % create multi-index
    muid = multiindex_create(multiindex_getPos(muid_offset), index_last);

    % compute sum `s`
    s = 0;
    isEnd = 0;
    while ~isEnd % loop over all multi-indices
        % get current cosition of multi-index
        index = multiindex_getPos(muid);

        % get 1-dimensional index
        index_1dim = multiindex_nDimTo1Dim_colMajor(index, ...
                                                    index_first, ...
                                                    index_length) + 1;

        % update sum `s`
        s = s + mask(index_1dim);

        % increment multi-index by two
        [muid,isEnd] = multiindex_increment_colMajor(muid, 2);
    end

    % scale if necessary
    if abs(s - scaleTo) > 10 * eps
        muid = multiindex_setPosToStart(muid);
        isEnd = 0;
        while ~isEnd % loop over all multi-indices
            % get current cosition of multi-index
            index = multiindex_getPos(muid);

            % get 1-dimensional index
            index_1dim = multiindex_nDimTo1Dim_colMajor(index, ...
                                                        index_first, ...
                                                        index_length) + 1;

            % norm mask coefficient
            mask(index_1dim) = scaleTo * mask(index_1dim) / s;

            % increment multi-index by two
            [muid,isEnd] = multiindex_increment_colMajor(muid, 2);
        end
    end

    % increment multi-index for offset
    [muid_offset,isEnd_offset] = multiindex_increment_rowMajor(muid_offset);
end

% end function
end
