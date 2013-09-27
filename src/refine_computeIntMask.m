%REFINE_COMPUTEINTMASK Compute mask coefficients of a refinable integral.
% Computes the mask coefficients of a refinable integral from the
% mask coefficients of the refinable functions in the integrand.
%
% Syntax:
% intMask = REFINE_COMPUTEINTMASK(dim, numFncs, mask, derivatives, ...
%                                 index_first, index_last, ...
%                                 IndexMask_first, IndexMask_last)
%
% Input:
% dim                dimension of the functions' domains
% numFncs            number of functions
% mask{}[]           cell-array with mask coefficients of the functions
% derivatives[]      derivatives of the functions
% index_first[]      beginning of multi-index of integral mask
% index_last[]       end of multi-index of integral mask
% IndexMask_first[]  beginning of multi-indices of all functions in integrand
% IndexMask_last[]   end of multi-indices of all functions in integrand
%
% Output:
% intMask            mask coefficients of a refinable integral
%
% See also: REFINE_COMPUTEINITIALVALUES, REFINE.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function intMask = refine_computeIntMask(dim, numFncs, mask, derivatives, ...
                                         index_first, index_last, ...
                                         IndexMask_first, IndexMask_last)

% set number of index elements
index_numElements = prod(abs(index_last - index_first) + 1);

% transpose multi-indices for mask coefficients of functions to simplify their
% usage in `refine_computeIntMaskEntry`
IndexMask_first = IndexMask_first';
IndexMask_last = IndexMask_last';

% create multi-index (input for function `refine_computeIntMaskEntry`)
muid_entry = multiindex_create(IndexMask_first(:,1), IndexMask_last(:,1));

% create multi-index for looping over mask coefficients of the integral
muid = multiindex_create(index_first, index_last);

% compute mask coefficients of the integral
intMask = zeros(1, index_numElements);
k = 1;
isEnd = 0;
while ~isEnd % loop over all multi-indices of `muid`
    % get current multi-index
    index = multiindex_getPos(muid);

    % compute one mask coefficient of integral
    intMask(k) = refine_computeIntMaskEntry(dim, numFncs, mask, index, ...
                                            muid_entry, ...
                                            IndexMask_first, IndexMask_last);

    % increment multi-index
    [muid,isEnd] = multiindex_increment_colMajor(muid);

    % increment `k`
    k = k + 1;
end

% end function
end




%REFINE_COMPUTEINTMASKENTRY Compute one entry of the integral mask coeff matrix.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function entry = refine_computeIntMaskEntry(dim, numFncs, mask, index, ...
                                            muid_entry, ...
                                            IndexMask_first, IndexMask_last)

% initiate return value that will be one mask coefficient of an integral
entry = 0;

% set current position of multi-index to first index
muid_entry = multiindex_setPosToFirst(muid_entry);

isEnd_entry = 0;
while ~isEnd_entry % loop over all multi-indices of first function in integral
    % get current position of multi-index
    index_entry = multiindex_getPos(muid_entry);

    % map multi-index `index_entry` to scalar index `index_entry_1dim`
    index_entry_1dim = multiindex_nDimTo1Dim_colMajor(index_entry, ...
                                                      IndexMask_first(:,1), ...
                                                      IndexMask_last(:,1)) + 1;

    % compute factor of a mask coeff. of the first function in the integral
    h = 1;
    for j = 2:numFncs % loop over all but first function constituted by `mask`
        subindex = index_entry - index((j-2)*dim+1:(j-1)*dim);
        if min(IndexMask_first(:,j) <= subindex) ...
           && min(subindex <= IndexMask_last(:,j))
            index_1dim = multiindex_nDimTo1Dim_colMajor( ...
                             subindex, ...
                             IndexMask_first(:,j), ...
                             IndexMask_last(:,j)) + 1;
            h = h * mask{j}(index_1dim);
        else
            h = 0;
            break
        end
    end

    % update `entry`
    entry = entry + mask{1}(index_entry_1dim) * h;

    % increment multi-index
    [muid_entry,isEnd_entry] = multiindex_increment_colMajor(muid_entry);
end

entry = entry / (2^dim);

% end function
end
