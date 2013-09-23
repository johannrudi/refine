%multiindex_isEnd  Check if multi-index reached its end.
%  Checks whether a multi-index structure has reached its end (if true then
%  current position of multi-index is beyond end).
%
%  Syntax:
%  isEnd = multiindex_isEnd(multiindex)
%
%  Input:
%  multiindex  multi-index structure as provided by `multiindex_create`
%
%  Output:
%  isEnd       `1` if multi-index reached its end, `0` otherwise
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

function isEnd = multiindex_isEnd(multiindex)

% set to `1` if current position is beyond end, `0` otherwise
isEnd = max(multiindex.pos > multiindex.index_last) && ...
        min(multiindex.pos >= multiindex.index_last);

% end function
end
