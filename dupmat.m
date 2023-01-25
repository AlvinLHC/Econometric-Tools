function d = dupmat(n)
% Returns Magnus and Neudecker's duplication matrix of size n
a = tril(ones(n));
i = find(a);
a(i) = 1:length(i);
a = a + tril(a,-1)';

j = a(:);

m = n*(n+1)/2;
d = zeros(n*n,m);
for r = 1:size(d,1)
    d(r, j(r)) = 1;
end