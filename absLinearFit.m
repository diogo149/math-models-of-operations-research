function [A] = absLinearFit(X, Y, ~)
%this function makes the tableau used to optimize
%absolute value functions with a linear fit
%
%we are trying to minimize the sum of abs(a*X+b-Y)
%
% let u_i - v_i = A X_i + b - Y_i
% thus u_i + v_i = |A X_i + b - Y_i|
%
% we also need two variables for A and b since negative values are allowed
% A = a1 - a2
% b = b1 - b2
%
% (this is the part of the code that would be changed if ever this was
%   not a linear fit)
%
% goal will be minimization of the sum of all u's and v's
%

%here we use A = A1 - w
%        and b = b1 - w
%to be more efficient
m = size(X,1);
%n = 1 (Y) + 1(A1) + 1(b1) + 1(w) + 
A = zeros(m+1,4+2*m);

A(2:end,1) = Y;
A(2:end,2) = X;
A(2:end,3) = ones(m,1);
A(2:end,4) = -X - ones(m,1);
A(1,5:end) = 1;
A(2:end, 4+(1:2:2*m)) = eye(m);
A(2:end, 4+(2:2:2*m)) = -eye(m);

%code to make it canonical form
if nargin == 3
    A = diag([1; 2*(Y >= 0)-1])*A;
    A(1,:) = A(1,:)-sum(A(2:end,:));
    e = 1;
    while e == 1
        [A i j e] = simplexRule(A);
    end
    exitFlag = 1;
    %deleting 0 rows
    z = all(~A,2);
    A(z,:) = [];
    x = effBFS(A);
end
            