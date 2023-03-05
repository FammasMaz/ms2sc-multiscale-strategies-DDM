% Algorithm 3 from Krylov Solvers PDF

% Initialization 
clear all; close all; clc;

addpath('base/');
run("data.m")

nbSub = truss.nbSub; % Number of Subdomains
nblocNodes = truss.nblocNodes; % Number of Elements per subdomain

% Primal Initialization
Sp_c = sparse(2*nbSub); % Initialization of Sp (concatenated)
bp_c = sparse(2*nbSub,1); % Initialization of bp (right hand side concatenated)
reshapeNodes = [1;]; 

% Dual Initialization
Sd_c = sparse(2*nbSub); % Initialization of Sp (concatenated)

% Loop for each sub-element 
% u_b (boundary) nodes
for i = 1:nbSub
    reshapeNodes = [reshapeNodes; i*nblocNodes+1];
end

u_b_c = sparse(length(reshapeNodes),1); % Unitialize the solution vector for u_b (0 for iteration 0) 
A_c = sparse(size(reshapeNodes,1), 2*nbSub); % Initialize A (Assembly Operator) - dimensions = Nb of subdomains +1 X Nb of subdomians*2
Ad_c = sparse(size(reshapeNodes,1), 2*nbSub); % Initialize concatenated dual (Assembly Operator) same dimension for A in a 1d element
ae = [1 0;0 1]; % First and last positions (map) - only works for 1D element
aed = [1 0;0 -1]; % First and last position (map) for dual - only works for 1D element
r_b = zeros(1,2*nbSub)'; % Residual 
R0_c  = zeros(1,nbSub)';
R_c = []; % Rigid Body modes concatenated

% Loop for each subdomain
for i=1:nbSub
    u_bs = u_b_c(i:i+1); % u_b on each subdomain
    [Sps, bps] = fem_k(truss, 0); % Gets the Sp for each subdomain
    Sds = pinv(full(Sps)); % Calculates the pseudoinverse of the Primal Shur - Equal the dual Shurr complement
    Rs = null(full(Sps),"r"); % Kernel of Primal Shure - Corresponds to rigid body modes
    R_c = blkdiag(R_c,Rs); % Rigid body modes concatenated
    if i == 1
        bps(1) = 0;
    end
    r_bs = bps - Sps*u_bs;
    R0s = sum(r_bs);
    R0_c(i) = R0s;
    r_b(2*i-1:2*i) = r_bs;
    
    % Concatenated assembly
    Sp_c(2*i-1:2*i,2*i-1:2*i) = Sps; % Stacks it on the general Sp (Stacked matrix)
    Sd_c(2*i-1:2*i,2*i-1:2*i) = Sds; % Stacks it on the concatenated Sb
    bp_c(2*i-1:2*i) = bps; % Gets the loading per element
    A_c(i:i+1,2*i-1:2*i) = A_c(i:i+1,2*i-1:2*i) + ae; % Constructs the Assembly Operator
    Ad_c(i:i+1,2*i-1:2*i) = Ad_c(i:i+1,2*i-1:2*i) + aed; % Constructs the Assembly dual Operator
end
d0 =  R0_c;

% Primal Approach 

% Assembled Quantities of primal
S = A_c*Sp_c*A_c'; % 
b = A_c*bp_c; %
r = A_c*r_b;

% Elimination of first rows and columns - First element has no displacement
Sr = S(2:end,2:end);
br = b(2:end);
u_b_0 = u_b_c(2:end); 
r_0 = r(2:end);

% Solution by Conjugate Gradient
[sol,it] = Conjugate_gradient(Sr,br,u_b_0,1e-5);

% Solution by Preconditioned Conjugate Gradient
A_tild = (A_c*A_c')\A_c;
P_Neumann = A_tild*Sp_c*A_tild';
P_Neumannr = P_Neumann(2:end,2:end);
P2 = diag(diag(Sr));

[sol2,it2,~] = gradPre(Sr,br,P_Neumannr,u_b_0,1e-7);

% Dual Approach 
bd_c = Sd_c*bp_c; % Dual right hand side

% Assembled Quantities
Sd = Ad_c*Sd_c*Ad_c';  
bd = Ad_c*bd_c; 
G  = Ad_c*R_c;
e = R_c'*bp_c;

Dual_Matrix = [Sd G;G' zeros(size(G,2))];
Dual_rhs  = [-bd; -e];
Dual_sol = Dual_Matrix\Dual_rhs;
Lambda = Dual_sol(1:size(Sd,2));
alpha = Dual_sol(size(Sd,2)+1:end);

% Interface displacement - currently not working
Lambda_c = Ad_c'*Lambda;
u_b_d = Sd_c*(bd_c + Lambda_c) + R_c*alpha;

clear all; close all; clc;
addpath('base/');

nbSub = 20; % number of subdomains
nbLocalElems = 10; % number of per domain elements

truss = mesher(nbSub, nbLocalElems); % generate truss fields

Sp_gen = 1 % primal problem


truss.DOF = (nbSub*(truss.nbNodes-1))+1;


reshapeNodes = [1;];
for i = 1:nbSub
    reshapeNodes = [reshapeNodes; i*truss.nblocNodes+1];
end

uord = sparse(length(reshapeNodes),1);

A = A_gen(reshapeNodes, truss.nbSub, Sp_gen);
[~, bp, Sp] = RS_gen(truss, truss.Fd, Sp_gen);

S = A*Sd*A'; % Overall
b = A*bp; % Overall


% Conjugate Gradient method
x0 = sparse(length(b), 1);
[Uord, iter, ~] = conjGradFunc(S, b, x0,200, 1e-5);
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
