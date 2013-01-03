function [A exitFlag] = effPhase0(A,tol)
%efficient phase 0 algorithm
if nargin == 1
    tol = 1e-8;
end
[m n] = size(A);
z = []; %storing the zero rows
for i = 2:m
    j = find(A(i,2:end)~= 0,1)+1;
    if size(j)
        A = piv(A,i,j,tol);
    else
        z = [z i];
    end
end

if ~size(z,1)
    exitFlag = 1;
elseif ~A(z,1)
    A(z,:) = [];
    exitFlag = 1;
else
    exitFlag = 0;
end