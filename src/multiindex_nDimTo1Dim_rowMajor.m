%multiindex_nDimTo1Dim_rowMajor  Map multi-index to scalar index (row-major).
%  Maps the current position of a multi-index to a unique one-dimensional
%  (scalar) value inside an interval beginning at `0`. This can be used to
%  address elements of a matrix in row-major order.
%
%  Syntax:
%  index_1dim = multiindex_nDimTo1Dim_rowMajor(index_ndim,
%                                              index_first,
%                                              index_last)
%
%  Input:
%  index_ndim[]   current position of a multi-index (integer array)
%  index_first[]  first index (integer array)
%  index_last[]   last index (integer array)
%
%  Output:
%  index_1dim     one-dimensional index (beginning at `0`)
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-06-25
%  ----------------------------------------------------------------------------

function index_1dim = multiindex_nDimTo1Dim_rowMajor(index_ndim, ...
                                                     index_first, ...
                                                     index_last)

% set length of multi-index
index_length = abs(index_last - index_first) + 1;

% set mapping from n-dim to 1-dim index
mapping = [index_length(2:end)' 1];
for k = length(index_ndim)-2:-1:1
    mapping(k) = mapping(k) * mapping(k+1);
end

% compute 1-dimensional index
index_1dim = mapping * (index_ndim - index_first);

% end function
end
