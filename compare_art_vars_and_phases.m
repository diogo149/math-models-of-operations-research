%scirpt used to compate artificial variable algorithm vs phase 0 and 1
%algorithm

%note: for accuracy measures, there may be mistakes with tableaus with
%multiple optimal solutions, but they are incredibly unlikely and thus not
%considered (potential improvement: use the optimal z instead of x)

function compare_art_vars_and_phases()
Ms = [5 10 20 50 100];
Ns_minus_Ms = [1 5 10 20 50 100 200];

time_matrix_art_vars = zeros(length(Ms),length(Ns_minus_Ms));
time_matrix_phases = zeros(length(Ms),length(Ns_minus_Ms));
accuracy_matrix_art_vars = zeros(length(Ms),length(Ns_minus_Ms));
accuracy_matrix_phases = zeros(length(Ms),length(Ns_minus_Ms));
for k1 = 1:length(Ms)
    m = Ms(k1);
    for k2 = 1:length(Ns_minus_Ms)
        n = Ns_minus_Ms(k2) + m;
        i = get_i(m,n);
        [sum_time_phases, sum_time_art_vars, accuracy_phases,...
            accuracy_art_vars]=compare(i,m,n);
        time_matrix_art_vars(k1,k2) = sum_time_art_vars;
        time_matrix_phases(k1,k2) = sum_time_phases;
        accuracy_matrix_art_vars(k1,k2) = accuracy_art_vars;
        accuracy_matrix_phases(k1,k2) = accuracy_phases;
    end
end
keyboard()
end

function i = get_i(m,n)
%this function is needed because as m and n get very large, we need the
%number of iterations to decrease

if m*n > 1000
    i = 10;
elseif m*n > 500
    i = 20;
elseif m*n > 100
    i = 100;
elseif m*n > 20
    i = 500;
else
    i = 1000;
end

end

function [sum_time_phases, sum_time_art_vars, ...
    accuracy_phases, accuracy_art_vars]=compare(i,m,n)
iterations = i;
tol = 1e-8;

time_art_vars = zeros(iterations,1);
time_phases = zeros(iterations,1);
accuracy_art_vars = 0;
accuracy_phases = 0;
for k = 1:iterations
    A = randomCF(m,n);
    sol = solve(A,tol);
    A = removeCF(A);
    if mod(k,2) == 0
        [sol_art_vars, t_art_vars] = time(A, @(x) effArtVars(x,tol));
        [sol_phases, t_phases] = time(A, @(x) effPhaseSimplex(x,tol));
    else
        [sol_phases, t_phases] = time(A, @(x) effPhaseSimplex(x,tol));
        [sol_art_vars, t_art_vars] = time(A, @(x) effArtVars(x,tol));
    end
    time_art_vars(k) = t_art_vars;
    time_phases(k) = t_phases;
    accuracy_art_vars=accuracy_art_vars+solutions_match(sol, sol_art_vars,tol);
    accuracy_phases=accuracy_phases+solutions_match(sol,sol_phases,tol);
end
sum_time_phases = sum(time_phases);
sum_time_art_vars = sum(time_art_vars);
end

function match = solutions_match(sol1, sol2,tol)
    match = 0;
    if all(size(sol1) == size(sol2))
        if max(abs(sol1-sol2)) < tol
            match = 1;
        end
    end
end

function [x t] = time(A, f)
    tic
    x = f(A);
    t = toc;
end

function x = solve(A,tol)
    e = 1;
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
        x = [];
    end
end