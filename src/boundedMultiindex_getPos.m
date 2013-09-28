%BOUNDEDMULTIINDEX_GETPOS Current position of bounded multi-index.
% Gets the current position of a bounded multi-index structure.
%
% Syntax:
% pos = BOUNDEDMULTIINDEX_GETPOS(multiindex)
%
% Input:
% multiindex  multi-index structure as provided by `boundedMultiindex_create`
%
% Output:
% pos[]       current position of multi-index (integer array)
%
% See also: BOUNDEDMULTIINDEX_CREATE, BOUNDEDMULTIINDEX_INCREMENT,
% BOUNDEDMULTIINDEX_ISEND
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function pos = boundedMultiindex_getPos(multiindex)

% get current position of bounded multi-index
pos = multiindex.pos;

% end function
end
