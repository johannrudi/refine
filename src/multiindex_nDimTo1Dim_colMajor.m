%MULTIINDEX_NDIMTO1DIM_COLMAJOR Map multi-index to scalar index (column-major).
% Maps the current position of a multi-index to a unique one-dimensional
% (scalar) value inside an interval beginning at `0`. This can be used to
% address elements of a matrix in column-major order.
%
% Syntax:
% index_1dim = MULTIINDEX_NDIMTO1DIM_COLMAJOR(index_ndim,
%                                             index_first,
%                                             index_last)
%
% Input:
% index_ndim[]   current position of a multi-index (integer array)
% index_first[]  first index (integer array)
% index_last[]   last index (integer array)
%
% Output:
% index_1dim     one-dimensional index (beginning at `0`)
%
% See also: MULTIINDEX_CREATE, MULTIINDEX_INCREMENT_COLMAJOR,
% MULTIINDEX_INCREMENT_ROWMAJOR, MULTIINDEX_NDIMTO1DIM_INCRBYPOWOFTWO,
% MULTIINDEX_NDIMTO1DIM_ROWMAJOR.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function index_1dim = multiindex_nDimTo1Dim_colMajor(index_ndim, ...
                                                     index_first, ...
                                                     index_last)

% set length of multi-index
index_length = abs(index_last - index_first) + 1;

% set mapping from n-dim to 1-dim index
mapping = [1 index_length(1:end-1)'];
for k = 3:length(index_ndim)
    mapping(k) = mapping(k) * mapping(k-1);
end

% compute 1-dimensional index
index_1dim = mapping * (index_ndim - index_first);

% end function
end
