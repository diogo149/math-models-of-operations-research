function [x A exitFlag] = effArtVars(A,tol)
%efficient artificial variables
if nargin == 1
    tol = 1e-8;
end
x = [];
exitFlag = 0;
[m n] = size(A);
%making b's positive
A = diag([1; 2*(A(2:end,1) >= 0)-1])*A;
%saving objective row
obj = A(1,:);
A(1,:) = -sum(A(2:end,:));
%adding identity columns
%here I skip the steps of pivoting on each artificial column
A = [A [zeros(1,m-1); eye(m-1)]];
%performing simplex algorithm
e = 1;
while e == 1
    [A i j e] = simplexRule(A,tol);
end
% if e^T y != 0, problem is infeasible
if abs(A(1,1)) > tol
    return
end

for j = n+1:n+m-1
    %only way that cost is 0 is if y is basic!
    if abs(A(1,j)) < tol
        %since we delete 0 rows later, we only need to pivot on nonzero row
        zCnt = 0;
        oneLoc = 0;
        for i = 2:m
            if abs(A(i,j)) < tol
                zCnt = zCnt+1;
            elseif abs(A(i,j)-1) < tol
                oneLoc = i;
            else
                break
            end
        end
        if oneLoc > 1 && zCnt == m-2
            for k = 2:n
                if A(oneLoc,k) > tol
                    A = piv(A,oneLoc,k,tol);
                    break
                end
            end
        end
    end
end
%deleting the Y rows
A(:,n+1:end) = [];
A(1,:) = obj;
%deleting 0 rows
z = all(~A,2);
z(1) = 0;
A(z,:) = [];
m = size(A,1);
%restoring the objective function
A(1,:) = obj;
%pivoting on the identity columns
for j = 2:n
    zCnt = 0;
    oneLoc = 0;
    for i = 2:m
        if abs(A(i,j)) < tol
            zCnt = zCnt+1;
        elseif abs(A(i,j)-1) < tol
            oneLoc = i;
        else
            break
        end
    end
    if oneLoc > 1 && zCnt == m-2
        A(1,:) = A(1,:) - A(1,j)*A(oneLoc,:);
    end
end
e = 1;
%pivot to optimality
while e == 1
    [A i j e] = simplexRule(A,tol);
end
if e == 2
    exitFlag = 1;
    %deleting 0 rows
    z = all(~A,2);
    z(1) = 0;
    A(z,:) = [];
    x = effBFS(A,tol);
else
    exitFlag = 2;
end