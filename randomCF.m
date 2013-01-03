function A = randomCF(m,n)
%generates random m+1 by n+1 canonical form matrix; have n > m
assert(n > m)
A = eye(m);
A = [zeros(1,m); A];
A = [A 2*rand(m+1,n-m)-1];
A = A(:,randperm(n));
A = [rand(m+1,1) A];