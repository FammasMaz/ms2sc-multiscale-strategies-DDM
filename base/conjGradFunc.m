function [sol, k] = conjGradFunc(A,b,x_0,tol);


    residus = - A * x_0 + b;
    p =  residus;
    k = 0;
    x = x_0;
    n = length(x_0);
    double_residus = residus'* residus;
    while norm(residus) > tol
        Ap  = A * p;
        alpha = double_residus / (p' * Ap);
        x = x + alpha * p;
        residus = residus - alpha * Ap;
        double_residus_act = residus'* residus;
        beta = double_residus_act / double_residus;
        double_residus = double_residus_act;
        p = residus + beta * p;
        k = k + 1;
    end
    sol = x;
end