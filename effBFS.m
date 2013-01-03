function [x S] = effBFS(A,tol)
%efficient get BFS
if nargin == 1
    tol = 1e-8;
end
x = [];
[m n] = size(A);
S = zeros(m-1,1);
for j = 2:n
    zCnt = 0;
    oneLoc = 0;
    for i = 1:m
        if abs(A(i,j)) < tol
            zCnt = zCnt+1;
        elseif abs(A(i,j)-1) < tol
            oneLoc = i;
        else
            break
        end
    end
    if oneLoc > 1 && zCnt == m-1
        S(oneLoc-1) = j-1;
    end
end
if min(S) == 0
    return
end
x = zeros(n-1,1);
x(S) = A(2:end,1);