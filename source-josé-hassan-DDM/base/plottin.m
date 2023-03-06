function plottin(truss, uxy)
    nElems = truss.nbElems;
    nNodes = truss.nbNodes;
    adjac  = zeros(truss.nbNodes); % adjacenecy matrix
    
    for i = 1:truss.nbElems
        ids = truss.elems(i, (1:end-1)); % Node IDs
        matPropind = truss.elems(i, end); % Material Prop
        x = ids(1); % Node 1
        y = ids(2); % Node 2
        adjac(x,y) = 1; adjac(y,x) = 1;
    end
 
    dXY = [truss.nodes 0.1*ones(size(uxy,1),1)] + [uxy 0.1*ones(size(uxy,1),1)]; % Coordinates of the deformed matrix
    
    
    % Plotting both the matrices
    figure
    gplot(adjac, [truss.nodes zeros(size(uxy,1),1)], "r");
    hold on;
    gplot(adjac, dXY, "b");
    ylim([-0.4 0.8]);
    title('Displacement');
    xlabel('Displacement (m)');
    legend('Original','Displaced');
    hold off;
end