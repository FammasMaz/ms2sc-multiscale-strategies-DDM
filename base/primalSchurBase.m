clear all; close all; clc;

% Notations
% Sps = subdomain S
% Sp = concatenated S
% S = Global S
% size of ub=S=b (ub has one less size because the first is 0 and we are
% using sparse)

% size of us =Sp=bp

addpath('base/');

run("data.m")
nbSub = truss.nbSub;
truss.DOF = (nbSub*(truss.nbNodes-1))+1;
nblocNodes = 10;
Sp = sparse(2*nbSub, 2*nbSub);
Sdp = sparse(2*nbSub, 2*nbSub);

bp = sparse(2*nbSub,1);

reshapeNodes = [1;];
for i = 1:nbSub
    reshapeNodes = [reshapeNodes; i*nblocNodes+1];
end
uord = sparse(length(reshapeNodes),1);

ae = [1 0; 0 1];
abare = []
R_c = [];
A = sparse(size(reshapeNodes,1), 2*(nbSub));
Fd = 10e5;
for i=1:nbSub
    [Sps, bps, Kii, Kib, fi] = fem_k(truss, sol)
    Spd = sparse(pinv(full(Sps)));

    rig = null(full(Sps),'r');

    if i==nbSub
        Sps = [0 0; 0 0];
        Spd = [0 0; 0 0];
        rig = [1];
    elseif i==1
        rig = null(5000);
    end
    R_c = blkdiag(R_c, rig);
    Sp(2*i-1:2*i,2*i-1:2*i) = Sps;
    Sdp(2*i-1:2*i,2*i-1:2*i) = Spd;
    bp(2*i-1:2*i) = bps;
    A(i:i+1,2*i-1:2*i) = A(i:i+1,2*i-1:2*i) + ae;    

end

A(:,1) = [];
A(1,:) = [];
A(end, :) = [];
A(:, end) = [];


bp(2*nbSub, 1) = Fd;


bcremOrd = zeros(length(bp), 1);
for n = 1:size(truss.BCD, 1)
  bcnode = truss.BCD(n,1);
  bcremOrd(bcnode) = truss.BCD(n,2);
end

Sp = Sp(~bcremOrd,~bcremOrd);
Sdp = Sdp(~bcremOrd,~bcremOrd);
Sp(1,1) = Sp(1,1)*4;
Sdp(1,1) = Sdp(1,1)*4;

bp = bp(~bcremOrd);
bd = Sp*bp;

b = A*bd;
G = A*R_c;
e = R_c'*bp;
S = A*Sp*A';
Sd = A*Sdp*A'
%%
lhd = [S G;G' zeros(size(G))];
rhd = [-b;-e];
sol = lhd\rhd;
lambda = sol(1:size(S,2));
alpha = sol(size(S,2)+1:end);
lambda_c = A'*lambda;
ub = Sp*(bd+lambda_c) + R_c*alpha;
ub = sparse(ub(2:2:end));


%%
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
[uii, u, uif, unf] = internalNodes_bak(truss, reshapeNodes, Sp, bp, uord)


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