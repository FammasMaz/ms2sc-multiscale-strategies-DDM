function [u, iter, residual] = conjGradPrimal(nbSub, nbLocElems, plt)
%   [u, iter, residual] = conjGradPrimal(nbSub, nbLocElems, plt)
% Solver for Primal Schur Problem using Conjugate Gradient Method
%   [u, uii, unf, S, hH] = primalSchur(nbSub, nbLocalElems, plt)
% Output :  u = total displacement vector
%           iter = number of iterations taken
%           residual = residual vector (each iteration)
% Input :   nbSub = number of subdomains
%           nbLocalElems = number of elements in subdomain
%           plt = to plot or not (0 or 1)

truss = mesher(nbSub, nbLocElems); % generate truss fields

Sp_gen = 1; % primal problem


truss.DOF = (nbSub*(truss.nbNodes-1))+1;


reshapeNodes = [1;];
for i = 1:nbSub
    reshapeNodes = [reshapeNodes; i*truss.nblocNodes+1];
end


A = A_gen(reshapeNodes, truss.nbSub, Sp_gen);
[~, bp, Sp] = RS_gen(truss, truss.Fd, Sp_gen);

S = A*Sp*A'; % Overall
b = A*bp; % Overall


% Conjugate Gradient method
x0 = sparse(length(b), 1);
[Uord, iter, residual] = conjGradFunc(S, b, x0,200, 1e-5);
ub_aug = sparse([0; Uord]);
[uii, u, uif, unf] = internalNodes(truss, truss.reshapeNodes, ub_aug);


% Rebuilding Truss
truss.nodes = [0:truss.h:truss.L*truss.nbSub]';
truss.nbNodes = length(truss.nodes);
truss.nbElems = truss.nbElems*truss.nbSub;
elems = [];
for i = 1:truss.nbElems
    elems = [elems; i i+1 1];
end
truss.elems = elems;

if plt == 1
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
title('Nodal Displacements calculated using Primal Schur (Conjugate Grad. Method)')
saveas(figure(1), fullfile('assets/nodal_u_primal_conjGrad.png'));
hold off;

plottin(truss, u)
title('Displaced Beam (Primal Schur: Conjugate Gradient)')% Plot without Reordering
saveas(figure(2), fullfile('assets/disp_beam_primal_conjGrad.png'));
end
