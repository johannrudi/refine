%boundedMultiindex_increment  Increment bounded multi-index (row-major).
%  Increments the current position of a bounded multi-index structure in
%  "row-major order," i.e., last entry of the multi-index gets incremented
%  first.
%
%  Syntax:
%  multiindex = boundedMultiindex_increment(multiindex)
%
%  Input:
%  multiindex  multi-index structure as provided by `boundedMultiindex_create`
%
%  Output:
%  multiindex  multi-index structure with incremented position
%
%  ----------------------------------------------------------------------------
%  Author:         Johann Rudi <johann@ices.utexas.edu>
%  ----------------------------------------------------------------------------

function multiindex = boundedMultiindex_increment(multiindex)

% set flag for incementation
incremented = 0;

% increment bounded multi-index (row-major)
for k = multiindex.dim:-1:1 % loop backward over each dimension
    if sum(abs(multiindex.pos - multiindex.index_first)) < multiindex.bound
                                       % if index magnitude is lower than bound
        % increment current multi-index position
        multiindex.pos(k) = multiindex.pos(k)+1;
        incremented = 1;
        break
    else
        % set k-th multi-index position on first value
        multiindex.pos(k) = multiindex.index_first(k);
    end
end

if ~incremented % if not incremented
    % set first entry of multi-index beyond bound
    multiindex.pos(1) = multiindex.index_first(1) + multiindex.bound + 1;
end

% end function
end
