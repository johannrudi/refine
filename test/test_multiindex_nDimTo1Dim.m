%TEST_MULTIINDEX_NDIMTO1DIM Test mappings from multi- to sclar indices.
% Tests the mappings of multi-indices to sclar (one-dimensional) indices by
% creating matrices and filling them with consecutive numbers in row-major
% order and column-major order.
%
% Syntax:
% TEST_MULTIINDEX_NDIMTO1DIM(l, m, n)
%
% Input:
% l, m, n  matrix dimensions
%
% ----------------------------------------------------------------------------
% Author:    Johann Rudi <johann@ices.utexas.edu>
% ----------------------------------------------------------------------------

function test_multiindex_nDimTo1Dim(l, m, n)

%% Test Mapping of 2-Dimensional Multi-Index

% create 2-dimensional matrices like (for m=3, n=4)
%   M2_colMajor =
%     1  4  7  10
%     2  5  8  11
%     3  6  9  12
%   M2_rowMajor =
%     1   2   3   4
%     5   6   7   8
%     9  10  11  12
M2_colMajor = zeros(m, n);
M2_rowMajor = zeros(n, m);
for k = 1:m*n
    M2_colMajor(k) = k;
    M2_rowMajor(k) = k;
end
M2_rowMajor = M2_rowMajor';

% create multi-index in order to access matrix rows and columns
index_first = [1 ; 1];
index_last = [m ; n];
muid = multiindex_create(index_first, index_last);

isEnd = 0;
while ~isEnd % loop over all rows and columns
    % get current multi-index
    index = multiindex_getPos(muid);

    % get 1-dimensional index
    index_1dim_colMajor = multiindex_nDimTo1Dim_colMajor(index, ...
                                                         index_first, ...
                                                         index_last) + 1;
    index_1dim_rowMajor = multiindex_nDimTo1Dim_rowMajor(index, ...
                                                         index_first, ...
                                                         index_last) + 1;

    % compare 1-dim index to matrix entry
    if M2_colMajor(index(1),index(2)) ~= index_1dim_colMajor || ...
       M2_rowMajor(index(1),index(2)) ~= index_1dim_rowMajor
        error('Test of mapping multi-index to 1-dim index failed.')
    end

    % increment multi-index
    [muid,isEnd] = multiindex_increment_rowMajor(muid);
end


%% Test Mapping of 3-Dimensional Multi-Index

% create 3-dimensional matrices like (for m=3, n=4)
M3_colMajor = zeros(l, m, n);
M3_rowMajor = zeros(n, m, l);
for k = 1:l*m*n
    M3_colMajor(k) = k;
    M3_rowMajor(k) = k;
end
%TODO finish setup of M3_rowMajor;

% create multi-index in order to access matrix rows and columns
index_first = [1 ; 1 ; 1];
index_last = [l ; m ; n];
muid = multiindex_create(index_first, index_last);

isEnd = 0;
while ~isEnd % loop over all rows and columns
    % get current multi-index
    index = multiindex_getPos(muid);

    % get 1-dimensional index
    index_1dim_colMajor = multiindex_nDimTo1Dim_colMajor(index, ...
                                                         index_first, ...
                                                         index_last) + 1;
    index_1dim_rowMajor = multiindex_nDimTo1Dim_rowMajor(index, ...
                                                         index_first, ...
                                                         index_last) + 1;

    % compare 1-dim index to matrix entry
    if M3_colMajor(index(1),index(2),index(3)) ~= index_1dim_colMajor
       %TODO compare matrix `M3_rowMajor`
        error('Test of mapping multi-index to 1-dim index failed.')
    end

    % increment multi-index
    [muid,isEnd] = multiindex_increment_rowMajor(muid);
end

% end function
end
