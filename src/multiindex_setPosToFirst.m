%MULTIINDEX_SETPOSTOFIRST Set multi-index to first.
% Reset current position of multi-index to first index.
%
% Syntax:
% multiindex = MULTIINDEX_SETPOSTOFIRST(multiindex)
%
% Input:
% multiindex  multi-index structure as provided by `multiindex_create`
%
% Output:
% multiindex  multi-index structure with new position
%
% See also: MULTIINDEX_CREATE, MULTIINDEX_GETPOS, MULTIINDEX_INCREMENT_COLMAJOR,
% MULTIINDEX_INCREMENT_ROWMAJOR, MULTIINDEX_ISEND.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function multiindex = multiindex_setPosToFirst(multiindex)

% reset position of multi-index
multiindex.pos = multiindex.index_first;

% end function
end
