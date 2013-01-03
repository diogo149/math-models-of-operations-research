function [P S]  = reinversion(A, S, origS, tol)
%performs reinversion
%M = P * A
%A = Q * M
%M = Q^-1 * A
%we want to find Q^-1 instead of P, which normally will have less error
%we do this by reconstructing our current basic sequence, taking at most
%m pivots
%
%when a loop happens in reinversion, two options:
%   1) change S and use new P
%   2) ignore reinversion and use old P
%trying option (1)

if nargin == 3
    tol = 1e-8;
end

m = size(A,1);
P = eye(m);

done = (S == origS) + 0.0;

k = 0;
oldDone = done;
while any(~done)
    k = k+1;
    if k > m - 1
        if all(oldDone == done)
            fprintf('loop in reinversion\n')
            S(~done) = origS(~done);
            return
        end
        k = 1;
        oldDone = done;
    end
    if done(k)
        continue
    end
    i = k+1;
    j = S(k)+1;
    p = P(i,:)*A(:,j);
    if abs(p) < tol
        continue        %not pivotable on zero
    end
    
    C = P*A(:,j);
    Q = eye(m);
    Q(:,i) = -C/p;
    Q(i,i) = 1/p;
    P = Q*P;
    
    done(k) = 1;
end