function [A] = piv(A,i,j,tol)
%minimalist version to pivot matrix A at the (i,j) element
assert(nargin > 2)
if nargin < 4
    tol = 1e-13;
end
p = A(i,j);
if abs(p) < tol
    A(i,j) = 0;
    return
end
I = eye(size(A,1));
I(:,i) = -A(:,j)/p;
I(i,i) = 1/p;
A = I*A;
A(abs(A) < tol) = 0;
A(i,j) = 1;