%BOUNDEDMULTIINDEX_CREATE Create bounded multi-index structure.
% Creates a structure array used for running loops over multi-dimensional
% indices that are bounded, i.e., the sum of the absolute values of each
% multi-index entry is bounded. The bounded multi-index structure has the
% following fields:
%
% multiindex.dim            dimension of the multi-index
% multiindex.index_first[]  first index (integer array of dimension `dim`)
% multiindex.bound          bound of multi-index
% multiindex.pos[]          current position of multi-index (integer array of
%                           dimension `dim`)
%
% Syntax:
% multiindex = BOUNDEDMULTIINDEX_CREATE(bound, index_first)
%
% Input:
% bound          upper bound for multi-index (greater or equal zero)
% index_first[]  first index (integer array)
%
% Output:
% multiindex     bounded multi-index structure
%
% See also: BOUNDEDMULTIINDEX_GETPOS, BOUNDEDMULTIINDEX_INCREMENT,
% BOUNDEDMULTIINDEX_ISEND
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function multiindex = boundedMultiindex_create(bound, index_first)

% set dimension
dim = length(index_first);

% create multi-index struct
multiindex = struct(...
    'dim', dim, ...
    'index_first', index_first(:), ...
    'bound', bound, ...
    'pos', index_first(:) ...
);

% end function
end
