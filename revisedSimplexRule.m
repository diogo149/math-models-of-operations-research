function [A e] = revisedSimplexRule(A,tol)
%performs repeated simplex rule with revised simplex method
%input:
%   A - canonical form simplex tableau
%output:
%   A - optimal form simplex tableau (if possible)
%   e - error flag: 0 for unbounded, 2 for optimal

if nargin == 1
    tol = 1e-8;
end

m = size(A,1);
P = eye(m);
[x origS] = effBFS(A,tol);
S = origS;
e = 1;
cnt = 0;
while e == 1
    cnt = cnt + 1;
    if cnt > m
        cnt = 0;
        [P S] = reinversion(A,S,origS,tol);
    end
    [P S e] = revisedSimplexStep(A,P,S);
end

A = P*A;