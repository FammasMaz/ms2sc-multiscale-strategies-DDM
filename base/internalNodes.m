function [ui, u, uif, unf] = internalNodes(truss, reshapeNodes, ub)

ui = [];
u = [];
unf = [];
uif = [];

X = truss.h*truss.nblocNodes*(truss.nbSub-2:truss.nbSub-1);
Y = full(ub(end-1:end));
yi = interp1(X, Y, truss.L*truss.nbSub, 'linear', 'extrap');
ub = [0;ub;yi];
ub
for i=1:truss.nbSub
    [~, ~, Kii, Kib, fi] = fem_k(truss, 0);
    uil = Kii\(fi - (Kib)*ub(i:i+1));
    ui = [ui; uil];
    u = [u; ub(i); Kii\(fi - (Kib)*ub(i:i+1));];
    unf = [unf; ub(i); zeros(length(uil),1)];
    uif = [uif; 0; Kii\(fi - (Kib)*ub(i:i+1));];
end
u
u = [u; ub(end)];
u
unf = [unf; ub(end)];
uif = [uif; 0];
end
