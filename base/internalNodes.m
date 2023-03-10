function [ui, u, uif, unf] = internalNodes(truss, reshapeNodes, ub)
    % Internal nodes and complete displacement vector GENERATOR

    %   [ui, u, uif, unf] = internalNodes(truss, reshapeNodes, ub)
    % !!! Only for uniform beam
    % Outputs:  Complete displacement vector `u`
    %          Internal nodal displacement `ui`
    %          Complete displacement vector with 0 as boundary nodes `uif`, useful for
    %               plots
    %          Boundary nodes displacement `unf`
    % Outputs: Kernel of the truss `R_c` representing the rigid body modes
    %               only works when Sp_gen == 0
    %          bp/bd (Sp_gen: 0 or 1 resp.) (concatenated) `bp`
    %          Sp/Sd (Sp_gen: 0 or 1 resp.) (concatenated) `Sp` 

    ui = [];
    u = [];
    unf = [];
    uif = [];

    for i=1:truss.nbSub-1
        [~, ~, Kii, Kib, fi] = fem_k(truss, 0, 0);
        uil = Kii\(fi - (Kib)*ub(i:i+1));
        ui = [ui; uil];
        u = [u; ub(i); Kii\(fi - (Kib)*ub(i:i+1));];
        unf = [unf; ub(i); zeros(length(uil),1)];
        uif = [uif; 0; Kii\(fi - (Kib)*ub(i:i+1));];
    end

    u = [u; ub(end)];


    % Last subdomain
    [~, ~, Kii, Kib, fi] = fem_k(truss, 0, 1);
    uie = Kii\(fi - (Kib)*ub(end));
    uife = [0; Kii\(fi - (Kib)*ub(end))];
    unf = [unf; ub(truss.nbSub)];
    u = [u; uie];
    uif = [uif; uife(1:end-1)];
end
