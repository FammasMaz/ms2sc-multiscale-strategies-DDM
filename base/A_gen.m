function A = A_gen(reshapeNodes, nbSub, Sp_gen)

if Sp_gen ==1
    ae = [1 0; 0 1];
else
    ae = [-1 0;0 1];
end
A = sparse(size(reshapeNodes,1), 2*(nbSub));
for i=1:nbSub
        A(i:i+1,2*i-1:2*i) = A(i:i+1,2*i-1:2*i) + ae;
end
A(:,1) = [];
A(1,:) = [];
A(end, :) = [];
A(:, end) = [];
end
