%multiindex_getPos  Current multi-index position.
%  Gets the current position of a multi-index structure.
%
%  Syntax:
%  pos = multiindex_getPos(multiindex)
%
%  Input:
%  multiindex  multi-index structure as provided by `multiindex_create`
%
%  Output:
%  pos[]       current position of multi-index (integer array)
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-01-03
%  ----------------------------------------------------------------------------

function pos = multiindex_getPos(multiindex)

% set current position of multi-index
pos = multiindex.pos;

% end function
end
