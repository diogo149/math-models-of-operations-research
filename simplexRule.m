function [A i j exitFlag] = simplexRule(A,tol)
%more efficient version of simplexStep: picks column with first negative c
%also picks first index with lowest b/A to prevent cycling
% this function does one step of the simplex rule on tableau A
% exitFlag = 2 if optimal form
% exitFlag = 1 if successful
% exitFlag = 0 if unbounded form
%this algorithm assumes A is in canonical form
if nargin == 1
    tol = 1e-8;
end
i=0;
j=0;
n = size(A,2);

for k = 2:n
    if A(1,k) < -tol
        j = k;
    end
end
if j == 0
    exitFlag = 2;
    return
end
if max(A(:,j))< tol
    exitFlag = 0;
    return
end

r = A(2:end, j);
r(r<0) = nan;
[m i] = min(A(2:end, 1)./r);
i = i(1)+1;
exitFlag = 1;
A = piv(A,i,j,tol);