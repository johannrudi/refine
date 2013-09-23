%boundedMultiindex_isEnd  Check if bounded multi-index reached its end.
%  Checks whether a bounded multi-index structure has reached its end (if true
%  then current position of multi-index is beyond bound).
%
%  Syntax:
%  isEnd = boundedMultiindex_isEnd(multiindex)
%
%  Input:
%  multiindex  multi-index structure as provided by `boundedMultiindex_create`
%
%  Output:
%  isEnd       `1` if multi-index reached its end, `0` otherwise
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  Last modified:  2012-06-22
%  ----------------------------------------------------------------------------

function isEnd = boundedMultiindex_isEnd(multiindex)

% set to `1` if current position is beyond bound, `0` otherwise
isEnd = sum(abs(multiindex.pos - multiindex.index_first)) > multiindex.bound;

% end function
end
