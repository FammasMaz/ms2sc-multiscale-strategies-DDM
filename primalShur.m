clear all; close all; clc;

addpath('base/');

run("data.m")
nbSub = 20;
truss.nbSub = nbSub;
nblocNodes = 10;
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
    Sp(2*i-1:2*i,2*i-1:2*i) = Sps;
    bp(2*i-1:2*i) = bps;
    us(2*i-1:2*i) = Sps\bps;
    A(i:i+1,2*i-1:2*i) = A(i:i+1,2*i-1:2*i) + ae;
end

bp(1,1) = 0;
S = A*Sp*A';
b = A*bp;

truss.BC = [1 1];
bcremOrd = zeros(length(b), 1);
for n = 1:size(truss.BC, 1)
  bcnode = truss.BC(n,1);
  bcremOrd(bcnode) = truss.BC(n,2);
end

rmKord = S(~bcremOrd,~bcremOrd); % Removing the corresponding rows and columns of bcrem
newFord = b(~bcremOrd); % Removing the corresponding rows of bcrem
% New U as Matrix Solution to [K]{u} = {F}
Uord = rmKord\newFord;
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
[uii, u] = internalNodes(truss, reshapeNodes, Sp, bp, uord)
plot(uxyOrd)
