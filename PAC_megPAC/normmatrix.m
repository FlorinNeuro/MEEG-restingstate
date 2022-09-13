function yi= normmatrix(xi)

% takes a single matrix and returns the matrices with rows normalized to a length of one.
cols = size(xi,2);
n = 1 ./ sqrt(sum(xi.*xi,2));
yi = xi .* n(:,ones(1,cols));
yi(~isfinite(yi)) = 1;