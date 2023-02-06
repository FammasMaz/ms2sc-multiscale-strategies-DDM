function [ui, u] = internalNodes(truss, reshapeNodes, Sp, bp, ub)

ui = [];
u = [];
for i=1:truss.nbSub
    [~, ~, Kii, Kib, fi] = fem_k(truss, 0);
    ui = [ui; Kii\(fi - (Kib)*ub(i:i+1))];
    u = [u; ub(i); Kii\(fi - (Kib)*ub(i:i+1));];
end
u = [u; ub(end)];
end

