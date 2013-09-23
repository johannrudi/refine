%multiindex_setPosToFirst  Set multi-index to first.
%  Reset current position of multi-index to first index.
%
%  Syntax:
%  multiindex = multiindex_setPosToFirst(multiindex)
%
%  Input:
%  multiindex  multi-index structure as provided by `multiindex_create`
%
%  Output:
%  multiindex  multi-index structure with new position
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-06-22
%  ----------------------------------------------------------------------------

function multiindex = multiindex_setPosToFirst(multiindex)

% reset position of multi-index
multiindex.pos = multiindex.index_first;

% end function
end
