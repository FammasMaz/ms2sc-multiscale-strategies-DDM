function [R_c, bp, Sp] = RS_gen(truss, Fd, Sp_gen)

% R_c is only correct when the Sp_gen ~=1
R_c = [];
Sp = sparse(2*truss.nbSub, 2*truss.nbSub);
bp = sparse(2*truss.nbSub,1);


for i=1:truss.nbSub
    [Sps, bps, ~, ~, ~] = fem_k_dual(truss, 0, Sp_gen);
    rig = null(full(Sps),'r');
    if i==truss.nbSub
        Sps = [0 0; 0 0];
        rig = [1];
    elseif i==1
        rig = null(5000);
    end
    R_c = blkdiag(R_c, rig);
    Sp(2*i-1:2*i,2*i-1:2*i) = Sps;
    bp(2*i-1:2*i) = bps;
end
% Inputting Force on the last node
bp(2*truss.nbSub, 1) = Fd;

% Removing the lines from S and b according to the number of number of
% boundary nodes
bcremOrd = zeros(length(bp), 1);
for n = 1:size(truss.BCD, 1)
  bcnode = truss.BCD(n,1);
  bcremOrd(bcnode) = truss.BCD(n,2);
end

Sp = Sp(~bcremOrd,~bcremOrd);
if Sp_gen ~= 1
    Sp(1,1) = Sp(1,1)*4;
end
bp = bp(~bcremOrd);
R_c = sparse(R_c);
end
