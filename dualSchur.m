function [u, uii, unf, S, hH] = dualSchur(nbSub, nbLocalElems, plt)
    % Solver for Dual Schur Problem using Direct Method
    %   [u, uii, unf, S, hH] = dualSchur(nbSub, nbLocalElems, plt)
    % Output :  u = total displacement vector
    %           uii = internal nodal displacment vector
    %           unf = boundary nodal displacement vector
    %           S = Sp assembled
    %           hH = h/H parameter
    % Input :   nbSub = number of subdomains
    %           nbLocalElems = number of elements in subdomain
    %           plt = to plot or not (0 or 1)

    % Generating Primal or Dual? Sp_gen = 1 for Primal, 0 for Dual
    Sp_gen = 0;
    % Notations
    % Sps = subdomain S
    % Sp = concatenated S
    % S = Global S
    % size of ub=S(:,1)=b (ub has one less size because the first is 0 and we are
    % using sparse)

    % size of us=Sp(:,1)=bp

    truss = mesher(nbSub, nbLocalElems);

    hH = truss.h/truss.L;
    % generate truss fields

    DOF = truss.DOF; % DOFs of the truss
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
    ub_aug = [0;ub];

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
    title('Nodal Displacements calculated using Dual Schur (Direct Method)')
    saveas(figure(1), fullfile('assets/nodal_u_dual_direct.png'));
    hold off;

    plottin(truss, u)
    title('Displaced Beam (Dual Schur: Direct)')% Plot without Reordering
    saveas(figure(2), fullfile('assets/disp_beam_dual_direct.png'));
    end
end
