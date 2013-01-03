function [P S e] = revisedSimplexStep(A,P,S,tol)
%revised simplex algorithm
%input:
%   A - original canonical form simplex tableau
%   P - pivot matrix so far
%   S - basic sequence
%output:
%   P - pivot matrix for next step
%   S - new basic sequence
%   e - error flag: 0 for unbounded, 1 for successful, 2 for optimal

if nargin == 3
    tol = 1e-8;
end

[m n] = size(A);
P1 = P(1,:); %first row of pivot matrix
optimal = 1;
for j = 2:n
    if ismember(j,S)
        continue
    end
    if P1 * A(:,j) < -tol
        optimal = 0;
        break
    end
end

if optimal
    e = 2;
    return
end

C = P*A(:,j); %pivot column in the new simplex tableau

loc = find(C>tol);
if isempty(loc)
    e = 0;
    return
end

B = P*A(:,1); %constant column
[temp i] = min(B(loc)./C(loc));
i = loc(i);

if temp < -tol
    fprintf('nonnegative B\n')
    keyboard()
end


p = C(i); %pivot element
Q = eye(m);
Q(:,i) = -C/p;
Q(i,i) = 1/p;
P = Q*P;

S(i-1) = j-1; %set new basic sequence
e = 1;