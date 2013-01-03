function [x A exitFlag] = effPhaseSimplex(A,tol)
%efficient algorithm using phase0 and phase1 algorithms
if nargin == 1
    tol = 1e-8;
end
exitFlag = 0;
x=[];
[A e] = effPhase0(A,tol);
if e == 0
    return
end
[A e] = effPhase1(A,tol);
if e == 0
    return
end
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
end