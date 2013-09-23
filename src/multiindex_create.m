%multiindex_create  Create multi-index structure.
%  Creates a structure array used for running loops over multi-dimensional
%  indices. The multi-index structure has the following fields:
%
%  multiindex.dim            dimension of the multi-index
%  multiindex.index_first[]  first index (integer array of dimension `dim`)
%  multiindex.index_last[]   last index (integer array of dimension `dim`)
%  multiindex.pos[]          current position of multi-index (integer/float
%                            array of dimension `dim`)
%
%  Syntax:
%  multiindex = multiindex_create(index_first, index_last)
%
%  Input:
%  index_first[]  first index (integer array)
%  index_last[]   last index (integer array of same dimension as `index_first`)
%
%  Output:
%  multiindex     multi-index structure
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

function multiindex = multiindex_create(index_first, index_last)

% set dimension
dim = length(index_first);

% check for consistency in dimensions of input indices
if dim ~= length(index_last)
    error('Index dimensions of first and last index are not equal.')
end

% create multi-index structure
multiindex = struct(...
    'dim', dim, ...
    'index_first', index_first(:), ...
    'index_last', index_last(:), ...
    'pos', index_first(:) ...
);

% end function
end
