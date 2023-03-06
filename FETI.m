function [u, it, Sdinv, Sp, Sd] = FETI(nbSub, nbLocalElems, iter, er_max, plt)
% Solver for Dual Schur Problem using FETI Method
%   [u, it, Spinv, Sp, Sd] = conjGradPreDual(nbSub, nbLocalElems, iter, er_max, plt)
% Output :  u = total displacement vector
%           it = number of iterations taken
%           Spinv = Spinv found from A and Sd_concat
%           Sp = assembled Sp
%           Sd = assembled Sd
% Input :   nbSub = number of subdomains
%           nbLocalElems = number of elements in subdomain
%           iter = maximum number of iterations
%           er_max = max error threshold
%           plt = to plot or not (0 or 1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables
truss = mesher(nbSub, nbLocalElems); % Generates truss mesh
Sp_gen = 0;
Fd = truss.Fd; % Force on the end node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assembling Stiffness Matrix Sp and Sd in concatenated form
[~, bp_concat, Sp_concat] = RS_gen(truss, Fd, Sp_gen);
[R, bd_test, Sd_concat] = RS_gen(truss, Fd, 0);

%% Assembling A and Abar matrices
A = A_gen(truss.reshapeNodes, truss.nbSub, Sp_gen);
Abar = A_gen(truss.reshapeNodes, truss.nbSub, 0);

%% Solving for ub
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

% Initial Conditions
u(:, 1) = -G*((G'*G)\e);
r(:, 1) = P'*bp;
z(:, 1) = Sdinv*r(:, 1);
d(:, 1) = z(:, 1);

it = 0; % Iteration counter
if norm(full(r)) > 10e-4 % Check if initial guess is good enough
    for i=1:iter
        % Check for convergence
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
        it = it+1;
    end
    ub = u(:, end);
else
    ub = u(:,1);
end

%% FETI Method
alpha_Fetti = Gb*(-bd-Sd*ub);
Lambda = Abar'*ub;
rig_mod = R*alpha_Fetti; % Rigidity Modulus
ub = Sd_concat*(bd_concat + Lambda) + rig_mod;

%% Reordering
ub = ub(2:2:end); 

%% Internal Nodes and plotting
ub_aug = [0;ub];
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

%Plotting and saving images
if plt == 1 % Plot with Reordering
figure
nonZeroUif = find(uif~=0);
nonZeroUnf = find(unf~=0);
plot(nonZeroUif,uif(nonZeroUif), 'go')
hold on;
plot(nonZeroUnf,unf(nonZeroUnf), 'rx')
legend('Internal Nodes','Boundary Nodes', 'Location', 'southeast');
xlabel('Node location on beam (m)')
ylabel('Node displacements (m)')
title('Nodal Displacements calculated using Dual Schur (FETTI Method)')
saveas(figure(1), fullfile('assets/nodal_u_dual_FETTI.png'));
hold off;

plottin(truss, u) 
title('Displaced Beam (Dual Schur: Fetti Method)')
saveas(figure(2), fullfile('assets/disp_beam_dual_FETTI.png'));

figure(3)
plot(rig_mod, 'rx');
title('Rigid Body Modes in Displacement');
xlabel('Rigid Body Index');
ylabel('Solution in Displacement');
saveas(figure(3), fullfile('assets/rig_beam_dual_FETTI.png'));
end