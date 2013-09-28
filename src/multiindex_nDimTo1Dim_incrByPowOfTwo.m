%MULTIINDEX_NDIMTO1DIM_INCRBYPOWOFTWO Map 2^k incr. multi-index to scalar.
% Maps the current position of a multi-index to a unique one-dimensional
% (scalar) value. The value that the multi-index is incremented by has to be
% a power of `2`.
%
% Syntax:
% index_1dim = MULTIINDEX_NDIMTO1DIM_INCRBYPOWOFTWO(index_ndim,
%                                                   index_first,
%                                                   index_last,
%                                                   incrementBy)
%
% Input:
% index_ndim[]   current position of a multi-index (integer/float array)
% index_first[]  first index (integer array)
% index_last[]   last index (integer array)
% incrementBy    scalar value that the multi-index is incremented by,
%                which has to be a power of `2`
%
% Output:
% index_1dim     one-dimensional index
%
% See also: MULTIINDEX_CREATE, MULTIINDEX_INCREMENT_COLMAJOR,
% MULTIINDEX_INCREMENT_ROWMAJOR, MULTIINDEX_NDIMTO1DIM_COLMAJOR,
% MULTIINDEX_NDIMTO1DIM_ROWMAJOR.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function index_1dim = multiindex_nDimTo1Dim_incrByPowOfTwo( ...
                          index_ndim, ...
                          index_first, ...
                          index_last, ...
                          incrementBy)

% Description of the algorithm: Let the dimension of the multi-index be `n = 3`,
% and
%   index_first = (a1, a2, a3)^T,
%   index_last  = (b1, b2, b3)^T,
%   index_ndim  = (c1, c2, c3)^T,
%   incrementBy = 2^k.
% We want to assign a unique scalar `m` to `index_ndim`. We compute
%   d1 = (|a1 - b1|/2^k + 1)^2,
%   d2 = (|a2 - b2|/2^k + 1)^1,
%   d3 = (|a3 - b3|/2^k + 1)^0,
% and then the result is
%   m = (d1, d2, d3)^T * (c1/2^k, c2/2^k, c3/2^k).

% compute one dimensional index
index_1dim = ( (abs(index_last - index_first) / incrementBy + 1)' ...
               .^ (length(index_ndim)-1:-1:0) ) ...
             * ((index_ndim - index_first) / incrementBy);

% end function
end
