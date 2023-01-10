function plottinOrd(uxyOrd)
    run('data.m');
    nElems = truss.nbElems;
    nNodes = truss.nbNodes;
    
    adjacOrd  = zeros(truss.nbNodes); % adjacenecy matrix
    
    for i = 1:truss.nbElems
        ids = truss.elems(i, (1:end-1)); % Node IDs
        idso = [find(truss.ordNodesids == ids(1)) find(truss.ordNodesids == ids(2))];
        matPropind = truss.elems(i, end); % Material Prop
        x = idso(1); % Node 1
        y = idso(2); % Node 2
        adjacOrd(x,y) = 1; adjacOrd(y,x) = 1;
    end
    
    dXYOrd = [truss.nodes(truss.ordNodesids) 0.2*ones(size(uxyOrd,1),1)] + [uxyOrd 0.2*ones(size(uxyOrd,1),1)]; % Coordinates of the deformed matrix
    
    % Plotting both the matrices
    figure
    gplot(adjacOrd, [truss.nodes zeros(size(uxyOrd,1),1)], "r--");
    hold on;
    gplot(adjacOrd, dXYOrd, "b");
    ylim([-0.4 0.8]);
    title('Displacement with Reordering');
    xlabel('Displacement (mm)');
    legend('Original','Displaced');
    hold off;
end