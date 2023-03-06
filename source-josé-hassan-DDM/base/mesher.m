function truss = mesher(nbSub, nblocNodes)% Variables
    % Generates a truss mesh, based on h and H
    % To generate one truss of uniform Hs (same number of elements in each
    % truss), input the nbSub and nblocNodes
    % To generate truss with non uniform Hs (same number of elements in each
    % truss), input the nbSub ==1 and nblocNodes, before putting this function
    % in a loop
    % Inputs = number of subdomains `nbSubs`
    %          number of smaller elemnts per subdomain `nblocNodes`
    % Output = truss mesh struct with various fields
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Number of subdomains

    truss.nbSub = nbSub;
    % Length
    L = 100/truss.nbSub; % Size of each sub domain
    truss.L = L;
    % Number of nbLocal
    truss.nblocNodes = nblocNodes;
    % Number of elements
    truss.nbElems = truss.nblocNodes;

    % Number of nodes
    truss.nbNodes = truss.nbElems+1;

    % Size of each Element
    h = L/truss.nbElems;
    truss.h = h;
    % Coordinates of the nodes (Position x, Position y)
    % Each row is a separate Node

    truss.nodes = [0:h:(truss.nbNodes-1)*h]';
    truss.nodesids = [1:truss.nbNodes]';

    % Dimension of Nodes
    truss.Dim = size(truss.nodes,2);

    % Connectivity matrix for elements (Node i, Node j, Material Prop Index)
    % Each row is a separate element
    elems = [];
    for i = 1:truss.nbElems
        elems = [elems; i i+1 1];
    end
    truss.elems = elems;


    % Boundary Conditions (Node id, x condition, y condition) 
    % 1 = fixed, 0 = free
    % Each row is a separate element
    truss.BC = [1 1];
    load_nodes = 2*truss.nbSub-1;
    truss.BCD = [1 1
                load_nodes 1];


    % Truss Loadings(Node id, Force in x, Force in y)
    truss.Fd = 10e5;
    truss.loads = [truss.nbNodes truss.Fd];


    % Truss Materials (Modulus, Surface area)
    % (Each row is a separate material)
    truss.mat = [2.1e11 10e-6];

    truss.reshapeNodes = [1;];
    for i = 1:truss.nbSub
        truss.reshapeNodes = [truss.reshapeNodes; i*truss.nblocNodes+1];
    end
    % DOF of truss
    truss.DOF = (truss.nbSub*(truss.nbNodes-1))+1;


    % Internal and ordered nodes ids and numbers used later in the code
    iinodes = setdiff(truss.nodesids,truss.reshapeNodes);
    truss.iinodes = length(iinodes);
    truss.ordNodesids = [iinodes; truss.reshapeNodes];

end