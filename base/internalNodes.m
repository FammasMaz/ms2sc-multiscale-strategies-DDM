function [ui, u, uif, unf] = internalNodes(truss, reshapeNodes, Sp, bp, ub)

ui = [];
u = [];
unf = [];
uif = [];
for i=1:truss.nbSub
    [~, ~, Kii, Kib, fi] = fem_k(truss, 0);
    uil = Kii\(fi - (Kib)*ub(i:i+1));
    ui = [ui; uil];
    u = [u; ub(i); Kii\(fi - (Kib)*ub(i:i+1));];
    unf = [unf; ub(i); zeros(length(uil),1)];
    uif = [uif; 0; Kii\(fi - (Kib)*ub(i:i+1));];
end
u = [u; ub(end)];
unf = [unf; ub(end)];
uif = [uif; 0];
end
