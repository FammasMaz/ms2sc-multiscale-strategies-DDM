function [ui, u, uif, unf] = internalNodes(truss, reshapeNodes, ub)

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
