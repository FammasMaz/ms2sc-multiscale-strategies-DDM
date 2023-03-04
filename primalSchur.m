clear all; close all; clc;

% Add base functions
addpath('base/');

% Generating Primal or Dual? Sp_gen = 1 for Primal, 0 for Dual
Sp_gen = 1;

%%%% Notations %%%%%%
% Sps = subdomain S
% Sp = concatenated S
% S = Global S
% size of ub=S=b (ub has one less size because the first is 0 and we are
% using sparse)
% size of us =Sp=bp
%%%%%%%%%%%%%%%%%%%%%

nbSub = 20; % number of subdomains
nbLocalElems = 10; % number of per domain elements

truss = mesher(nbSub, nbLocalElems); % generate truss fields

DOF = truss.DOF; % DOFs of the truss

% Generate A and A_bar assembly operators
A = A_gen(truss.reshapeNodes, truss.nbSub, 1); 
Abar = A_gen(truss.reshapeNodes, truss.nbSub, Sp_gen);

% Generate bp and Sp, kernel not needed
[~, bp, Sp] = RS_gen(truss, truss.Fd, Sp_gen);


% Sp Assembled
S = A*Sp*A';
b = A*bp;

% U using inverse of S
ub = S\b;
ub_aug = [0;ub];

% Internal Nodes calculations
[uii, u, uif, unf] = internalNodes(truss, truss.reshapeNodes, ub_aug);

% Reassigning some truss fields
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
legend('Internal Nodes','Boundary Nodes', 'Location', 'southeast');
xlabel('Node location on beam (m)')
ylabel('Node displacements (m)')
title('Nodal Displacements calculated using Primal Schur (Direct Method)')
saveas(figure(1), fullfile('assets/nodal_u_primal_direct.png'));
hold off;

plottin(truss, u)
title('Displaced Beam (Primal Schur: Direct)')% Plot without Reordering
saveas(figure(2), fullfile('assets/disp_beam_primal_direct.png'));
