clear all; close all; clc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables
addpath('base/');
run('data.m'); % Generates truss mesh
Sp_gen = 1;
iter = 0;
Fd = 10e5; % Force on the end node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




[~, bp_concat, Sp_concat] = RS_gen(truss, Fd, Sp_gen);
[R, bd_test, Sd_concat] = RS_gen(truss, Fd, 0);


A = A_gen(truss.reshapeNodes, truss.nbSub, Sp_gen);
Abar = A_gen(truss.reshapeNodes, truss.nbSub, 0);

Atil = (Abar*Abar')\Abar;

G = (Abar*R);

Sp = A*Sp_concat*A';
Sd = Abar*Sd_concat*Abar';
bp = A*bp_concat;
bd_concat = Sd_concat*bp_concat;
bd = Abar*bd_concat;
e = R'*bp_concat;


Sdinv = Atil*Sp_concat*Atil';

Gb = ((G'*G)\G');
P = sparse(eye(size(Sp, 1)))  - (G*Gb);

% Initialize Sparse
u = sparse(size(Sp, 1), iter+1);
p = sparse(size(Sp, 1), iter+1);
z = sparse(size(Sp, 1), iter+1);
d = sparse(size(Sp, 1), iter+1);
r = sparse(size(Sp, 1), iter+1);
alpha = sparse(size(Sp, 1), iter+1);



u(:, 1) = -G*((G'*G)\e);
r(:, 1) = P'*bp;
z(:, 1) = Sdinv*r(:, 1);
d(:, 1) = z(:, 1);


if norm(full(r)) > 10e-4

    for i=1:iter
        
        p(:, i) = P'*Sp*d(:, i); 
        alpha(:, i) = r(:, i)'*d(:, i)/(d(:, i)'*p(:, i));
        u(:, i+1) = u(:, i) + alpha(:, i).*d(:,i);
        r(:,i+1) = r(:,i) - alpha(:, i).*p(:,i);
        z(:,i+1) = Sdinv*r(:,i+1);
        beta = 0;
        for j = 1:i
            beta = beta - z(:,i+1)'*p(:,j)/(d(:,j)'*p(:,j));
        end
        d(:,i+1) = z(:,i+1) + beta.*d(:,i);
        
    end
    ub = u(:, end);
else
    ub = u(:,1);
end

alpha_Fetti = Gb*(-bd-Sd*ub);
Lambda = Abar'*ub;
u_b = Sd_concat*(bd_concat + Lambda) + R*alpha_Fetti;