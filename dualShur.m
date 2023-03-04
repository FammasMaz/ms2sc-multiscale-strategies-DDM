clear all; close all; clc;

% Generating Primal or Dual? Sp_gen = 1 for Primal, 0 for Dual
Sp_gen = 0;
% Notations
% Sps = subdomain S
% Sp = concatenated S
% S = Global S
% size of ub=S=b (ub has one less size because the first is 0 and we are
% using sparse)

% size of us =Sp=bp

addpath('base/');

run("data.m");
nbSub = truss.nbSub;
truss.DOF = (truss.nbSub*(truss.nbNodes-1))+1;
nblocNodes = truss.nblocNodes;


Fd = 10e5;

% Generate A and A_bar
A = A_gen(truss.reshapeNodes, truss.nbSub, 1);
Abar = A_gen(truss.reshapeNodes, truss.nbSub, Sp_gen);

[R_c, bp, Sp] = RS_gen(truss, Fd, Sp_gen);


% b dual
bd = Sp*bp;

% Assembled b
b = Abar*bd;

G = Abar*R_c;

e = R_c'*bp;

% Assembled Sp
S = Abar*Sp*Abar';
lhd = [S G;G' zeros(size(G))];
rhd = [-b;-e];
sol = lhd\rhd;
lambda = sol(1:size(S,2));
alpha = sol(size(S,2)+1:end);
lambda_c = Abar'*lambda;
ub = Sp*(bd+lambda_c) + R_c*alpha;
ub = ub(2:2:end); % Removing the repeating elements


[uii, u, uif, unf] = internalNodes(truss, truss.reshapeNodes, ub);

% Rebuilding Truss
truss.nodes = [0:truss.h:truss.L*truss.nbSub]';
truss.nbNodes = length(truss.nodes);
truss.nbElems = truss.nbElems*truss.nbSub;
elems = [];
for i = 1:truss.nbElems
    elems = [elems; i i+1 1];
end
truss.elems = elems;


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