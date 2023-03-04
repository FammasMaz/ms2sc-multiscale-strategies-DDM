clear all; close all; clc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables
run('data.m'); % Generates truss mesh
Sp_gen = 1;
iter = 3;
Fd = 10e5; % Force on the end node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[~, bp_concat, Sp_concat] = RS_gen(truss, Fd, Sp_gen);
[R, ~, Sd_concat] = RS_gen(truss, Fd, 0);


A = A_gen(truss.reshapeNodes, truss.nbSub, Sp_gen);
Abar = A_gen(truss.reshapeNodes, truss.nbSub, 0);

Atil = (A*A')\A;

G = (Atil*R);

Sp = A*Sp_concat*A';
bp = A*bp_concat;

Spinv = Atil*Sd_concat*Atil';

Gb = ((G'*Sp*G)\G');
P = sparse(eye(size(Sp, 1)))  - (G*Gb*Sp);

% Initialize Sparse
u = sparse(size(Sp, 1), iter+1);
p = sparse(size(Sp, 1), iter+1);
z = sparse(size(Sp, 1), iter+1);
d = sparse(size(Sp, 1), iter+1);
r = sparse(size(Sp, 1), iter+1);
alpha = sparse(size(Sp, 1), iter+1);


u(:, 1) = G*Gb*bp;
r(:, 1) = P'*bp;
z(:, 1) = Spinv*r(:, 1);
d(:, 1) = z(:, 1);
if norm(full(r))> 10e-4
    for i=1:iter
        p(:, i) = P'*Sp*d(:, i); 
        alpha(:, i) = r(:, i)'*d(:, i)/(d(:, i)'*p(:, i));
        u(:, i+1) = u(:, i) + alpha(:, i).*d(:,i);
        r(:,i+1) = r(:,i) - alpha(:, i).*p(:,i);
        z(:,i+1) = Spinv*r(:,i+1);
        beta = 0;
        for j = 1:i
            beta = beta - z(:,i+1)'*p(:,j)/(d(:,j)'*p(:,j));
        end
        d(:,i+1) = z(:,i+1) + beta.*d(:,i);
    end
else
    ub = u(:,1);
end









