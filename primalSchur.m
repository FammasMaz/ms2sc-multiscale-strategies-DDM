clear all; close all; clc;

% Generating Primal or Dual? Sp_gen = 1 for Primal, 0 for Dual
Sp_gen = 1;
% Notations
% Sps = subdomain S
% Sp = concatenated S
% S = Global S
% size of ub=S=b (ub has one less size because the first is 0 and we are
% using sparse)

% size of us =Sp=bp

addpath('base/');

run('data.m');
truss.DOF = (truss.nbSub*(truss.nbNodes-1))+1;

Fd = 10e5;

% Generate A and A_bar
A = A_gen(truss.reshapeNodes, truss.nbSub, 1);
Abar = A_gen(truss.reshapeNodes, truss.nbSub, Sp_gen);

[R_c, bp, Sp] = RS_gen(truss, Fd, Sp_gen);


% Sp Assembled
S = A*Sp*A';
b = A*bp;

ub = S\b;
ub_aug = [0;ub;47];
[uii, u, uif, unf] = internalNodes(truss, truss.reshapeNodes, ub_aug)
%Plotting
figure
nonZeroUif = find(uif~=0);
nonZeroUnf = find(unf~=0);
plot(nonZeroUif,uif(nonZeroUif), 'go')
hold on;
plot(nonZeroUnf,unf(nonZeroUnf), 'rx')
legend('Internal Nodes','Boundary Nodes');
hold off;
plottin(truss, u)