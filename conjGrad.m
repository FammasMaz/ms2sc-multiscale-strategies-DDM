clear all; close all; clc;

addpath('base/');

run("data.m")
nbSub = truss.nbSub;
truss.DOF = (nbSub*(truss.nbNodes-1))+1;
nblocNodes = truss.nblocNodes;
Sp = sparse(2*nbSub);
bp = sparse(2*nbSub,1);

reshapeNodes = [1;];
for i = 1:nbSub
    reshapeNodes = [reshapeNodes; i*nblocNodes+1];
end

uord = sparse(length(reshapeNodes),1);
A = sparse(size(reshapeNodes,1), 2*nbSub);
ae = [1 0;0 1];
for i=1:nbSub
    [Sps, bps, ~, ~, ~] = fem_k(truss, 0);
    if i==1
        bps(1) = 0;
    end
    Sp(2*i-1:2*i,2*i-1:2*i) = Sps; % Subdomain level
    bp(2*i-1:2*i) = bps; % Subdomain level
    A(i:i+1,2*i-1:2*i) = A(i:i+1,2*i-1:2*i) + ae;
end
S = A*Sp*A'; % Overall
b = A*bp; % Overall

truss.BC = [1 1];
bcremOrd = zeros(length(b), 1);
for n = 1:size(truss.BC, 1)
  bcnode = truss.BC(n,1);
  bcremOrd(bcnode) = truss.BC(n,2);
end

rmKord = S(~bcremOrd,~bcremOrd); % Removing the corresponding rows and columns of bcrem
newFord = b(~bcremOrd); % Removing the corresponding rows of bcrem
% New U as Matrix Solution to [K]{u} = {F}
%% Conjugate Gradient method
x0 = sparse(length(newFord),1);
[Uord, iter, ~] = conjGradFunc(rmKord, newFord, x0,200, 1e-5);

ub = Uord;
% Rentering the previosuly removed nodal data, ordered
j = 1;
for i = 1:size(uord, 1)
  if bcremOrd(i) == 1
    uord(i) = 0; % 0 strain because of fixed support
  else
    uord(i) = Uord(j); % value from U as no fixed support
    j = j+1;
  end
end
for i=1:length(b)
    uxyOrd(i,1) = uord(i);
end
[uii, u, uif, unf] = internalNodes(truss, reshapeNodes, Sp, bp, uord)


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