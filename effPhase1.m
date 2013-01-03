function [A exitFlag] = effPhase1(A,tol)
%efficient phase 1 algorithm
if nargin == 1
    tol = 1e-8;
end
exitFlag = 0;
n = size(A,2);

posB = find(A(2:end,1) >= 0); % vector of nonnegative b's
posB = posB + 1;
negB = find(A(2:end,1) < 0); % vector of negative b's
negB = negB + 1;

if ~any(posB)     %when there are no positive B's
    j = 2;
    i = negB(1);
    while j <= n
        if A(i,j) < 0
            break
        end
        j = j+1;
    end
    if j > n
        return  %infeasible form, negative b with nonnegativs a's
    end
    A = piv(A,i,j,tol);
    posB = negB(1);
    negB(1) = [];
end
while negB
    focus = negB(1);
    %checkpoint
    e = 1;
    while e
        [A e] = subproblemSimplex(A, focus, posB, tol);
    end
    if A(focus,1) < 0
        return
    end
    posB = [posB; focus];
    negB(1) = [];
end
%delete blank rows here?
exitFlag = 1;
end

function [A exitFlag] = subproblemSimplex(A,focus,Bs,tol)
%exitFlag = 0 means done
%exitFlag = 1 means not done
if A(focus,1) >= 0
    exitFlag = 0;
    return
end
j=0;
m = size(Bs,1);
n = size(A,2);
for k = n:-1:2 %change this to go through Bs!!!
    if A(focus,k) < 0
        j = k; %have to get lowest j instead of last j!
        if max(A(Bs,j))<= 0
            %unbounded supproblem! pivot on first negative
            A = piv(A,focus,j,tol);
            exitFlag = 0;
            return
        end
    end
end
if j == 0 %if optimal form, return
    exitFlag = 0;
    return
end

r = A(Bs, j);
r(r<0) = nan;
[m i] = min(A(Bs, 1)./r);
i = Bs(i);
exitFlag = 1;
A = piv(A,i,j,tol);
end