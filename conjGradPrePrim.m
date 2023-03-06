function [u, it, Spinv, Sp, Sd] = conjGradPrePrim(nbSub, nbLocalElems, iter, er_max, plt)
    % Solver for Primal Schur Problem using Preconditioned Conjugate Gradient Method
    %   [u, it, Spinv, Sp, Sd] = conjGradPrePrim(nbSub, nbLocalElems, iter, er_max, plt)
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
    Sp_gen = 1;
    Fd = truss.Fd; % Force on the end node
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %% Assembling Stiffness Matrix
    [~, bp_concat, Sp_concat] = RS_gen(truss, Fd, Sp_gen);
    [R, ~, Sd_concat] = RS_gen(truss, Fd, 0);

    %% Generation of A and Abar
    A = A_gen(truss.reshapeNodes, truss.nbSub, Sp_gen);
    Abar = A_gen(truss.reshapeNodes, truss.nbSub, 0);

    %% Solving for ub
    Atil = (A*A')\A;
    G = (Atil*R);
    Sp = A*Sp_concat*A';
    bp = A*bp_concat;
    Spinv = Atil*Sd_concat*Atil';

    Gb = ((G'*Sp*G)\G');
    P = sparse(eye(size(Sp, 1)))  - (G*Gb*Sp);

    % Initialize Sparse
    u = sparse(size(Sp, 1), iter+1);
    p = sparse(size(Sp, 1), iter+1);
    z = sparse(size(Sp, 1), iter+1);
    d = sparse(size(Sp, 1), iter+1);
    r = sparse(size(Sp, 1), iter+1);
    alpha = sparse(size(Sp, 1), iter+1);

    % Initial Guess
    u(:, 1) = G*Gb*bp;
    r(:, 1) = P'*bp;
    z(:, 1) = Spinv*r(:, 1);
    d(:, 1) = z(:, 1);

    it = 0; % Iteration counter
    if norm(full(r))> er_max % Check if initial guess is good enough
        for i=1:iter
            p(:, i) = P'*Sp*d(:, i); 
            alpha(:, i) = r(:, i)'*d(:, i)/(d(:, i)'*p(:, i));
            u(:, i+1) = u(:, i) + alpha(:, i).*d(:,i);
            r(:,i+1) = r(:,i) - alpha(:, i).*p(:,i);
            z(:,i+1) = Spinv*r(:,i+1);
            beta = 0;
            for j = 1:i
                beta = beta - z(:,i+1)'*p(:,j)/(d(:,j)'*p(:,j));
            end
            d(:,i+1) = z(:,i+1) + beta.*d(:,i);
            it = it +1;
        end
        ub = u(:,end);
    else
        ub = u(:,1);
    end

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

    % Plotting 
    if plt == 1
    figure
    nonZeroUif = find(uif~=0);
    nonZeroUnf = find(unf~=0);
    plot(nonZeroUif,uif(nonZeroUif), 'go')
    hold on;
    plot(nonZeroUnf,unf(nonZeroUnf), 'rx')
    legend('Internal Nodes','Boundary Nodes', 'Location', 'southeast');
    xlabel('Node location on beam (m)')
    ylabel('Node displacements (m)')
    title('Nodal Displacements calculated using Primal Schur (Preconditioned Conjugate Grad.)')
    saveas(figure(1), fullfile('assets/nodal_u_primal_conjGradPre.png'));
    hold off;

    plottin(truss, u)
    title('Displaced Beam (Primal Schur: Pre. Conjugate Gradient)')% Plot without Reordering
    saveas(figure(2), fullfile('assets/disp_beam_primal_conjGradPre.png'));
    end
end