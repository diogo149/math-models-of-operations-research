function [A] = remove0Rows(A)
%removes empty rows of tableau
%note: over half the time consumed from this function is used calling it
%so you're better off copying it if you want maximum performance

%deleting 0 rows
z = all(~A,2);
z(1) = 0;
A(z,:) = [];

%make sure to recalculate size afterwards if needed