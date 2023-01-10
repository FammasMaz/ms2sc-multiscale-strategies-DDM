%% Dataset

% Variables
L = 100/20;

% Number of subdomains
truss.nbSub = 1;
% Number of nbLocal
truss.nblocNodes = 10;
% Number of elements
truss.nbElems = truss.nbSub*truss.nblocNodes;

% Number of nodes
truss.nbNodes = truss.nbElems+1;

% Size of each Element
h = L/truss.nbElems;

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

% Truss Loadings(Node id, Force in x, Force in y)
truss.loads = [truss.nbNodes 1000000];


% Truss Materials (Modulus, Surface area)
% (Each row is a separate material)
truss.mat = [2.1e11 10e-6];

truss.reshapeNodes = [1;]
for i = 1:truss.nbSub
    truss.reshapeNodes = [truss.reshapeNodes; i*truss.nblocNodes+1]
end

iinodes = setdiff(truss.nodesids,truss.reshapeNodes);
truss.iinodes = length(iinodes);
truss.ordNodesids = [iinodes; truss.reshapeNodes];
