function [x0, er, i] = conjGradPreFunc(A, b, m,x0, er_max, P)
% Conjugate Gradient Preconditioned algorithm

i = 0;
r0 = -A*x0 + b;
z0 = P\r0;
p0 = z0;
error = norm(r0)/norm(b);
er = [];
while i<m && error>er_max
    alpha = (r0'*z0)/(p0'*A*p0);
    x0 = x0 + alpha*p0;
    rnew = r0 - alpha*(A*p0);
    znew = P\rnew;
    beta = (rnew'*znew)/(r0'*z0);
    p0 = znew + beta*p0;
    error = norm(rnew)/norm(b);
    er = [er error];
    i = i + 1;
    z0 = znew;
    r0 = rnew;
end
end
