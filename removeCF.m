function A = removeCF(A)
%takes in a canonical form tableau, and converts it to noncanonical form
%(used for testing)
[m n] = size(A);
A = [ [1;zeros(m-1,1)] rand(m,m-1) ]*A;