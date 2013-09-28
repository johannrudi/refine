%MULTIINDEX_GETPOS Current multi-index position.
% Gets the current position of a multi-index structure.
%
% Syntax:
% pos = MULTIINDEX_GETPOS(multiindex)
%
% Input:
% multiindex  multi-index structure as provided by `multiindex_create`
%
% Output:
% pos[]       current position of multi-index (integer array)
%
% See also: MULTIINDEX_CREATE, MULTIINDEX_INCREMENT_COLMAJOR,
% MULTIINDEX_INCREMENT_ROWMAJOR, MULTIINDEX_ISEND, MULTIINDEX_SETPOSTOFIRST.
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function pos = multiindex_getPos(multiindex)

% set current position of multi-index
pos = multiindex.pos;

% end function
end
