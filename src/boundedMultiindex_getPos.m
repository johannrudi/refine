%boundedMultiindex_getPos  Current position of bounded multi-index.
%  Gets the current position of a bounded multi-index structure.
%
%  Syntax:
%  pos = boundedMultiindex_getPos(multiindex)
%
%  Input:
%  multiindex  multi-index structure as provided by `boundedMultiindex_create`
%
%  Output:
%  pos[]       current position of multi-index (integer array)
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-01-03
%  ----------------------------------------------------------------------------

function pos = boundedMultiindex_getPos(multiindex)

% get current position of bounded multi-index
pos = multiindex.pos;

% end function
end
