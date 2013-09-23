%multiindex_increment_rowMajor  Increment multi-index (row-major).
%  Increments the current position of a multi-index structure in
%  "row-major order," i.e., last entry of the multi-index gets incremented
%  first.
%
%  Syntax:
%  multiindex = multiindex_increment_rowMajor(multiindex, [incrementBy])
%  [multiindex,isEnd] = multiindex_increment_rowMajor(multiindex, [incrementBy])
%
%  Input:
%  multiindex   multi-index structure as provided by `multiindex_create`
%  incrementBy  scalar value between zero and one that the multi-index is
%               incremented by (optional input, default is `1`)
%
%  Output:
%  multiindex   multi-index structure with incremented position
%  isEnd        `1` if multi-index reached its end, i.e., multi-index position
%               is beyond end, `0` otherwise (optional output)
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

function [multiindex,varargout] = multiindex_increment_rowMajor(multiindex, ...
                                                                incrementBy)

%%  Assume: Multi-Index Position Is Not Beyond End Index (Check Follows Below)

% check optional input
if nargin < 2
    incrementBy = 1;
end

% increment multi-index (row-major)
for k = multiindex.dim:-1:1 % loop backward over each dimension
    if (multiindex.pos(k) + incrementBy) <= multiindex.index_last(k)
                                     % if not beyond last index after next step
        % increment current multi-index position
        multiindex.pos(k) = multiindex.pos(k) + incrementBy;

        % set output value that multi-index position is not beyond end index
        if nargout == 2
            varargout = {0};
        end

        % return and terminate function
        return
    else
        % set k-th multi-index position on first value
        multiindex.pos(k) = multiindex.index_first(k);
    end
end


%% If Multi-Index Position Was Beyond End Index

% set multi-index position to end
multiindex.pos = multiindex.index_last;

% set multi-index position beyond end
multiindex.pos(end) = multiindex.pos(end) + incrementBy;

% set output value that multi-index position is beyond end index
if nargout == 2
    varargout = {1};
end

% end function
end
