function [Sps, bps, Kii, Kib, fi] = fem_k(truss, sol, last)
%% Add Paths

addpath("base/");
%% Initialize Dataset

%% Initialize Matrices

DoF = truss.Dim*truss.nbNodes;
Kord =sparse(DoF, DoF);
length = zeros(truss.nbElems,1);
uord = sparse(DoF,1);
Ford = sparse(DoF,1);

%% Local Calculations

ordNodesids = truss.ordNodesids;
if last == 1
    pop = find(ordNodesids == truss.nblocNodes+1);
    ordNodesids(pop) = [];
    ins = find(ordNodesids == truss.nblocNodes);
    ordNodesids = [ordNodesids(1:ins); truss.nblocNodes+1; ordNodesids(ins+1:end)];
end

for i=1:truss.nbElems
    ids = truss.elems(i, (1:end-1)); % Node IDs
    idsOrd = [find(ordNodesids == ids(1)) find(ordNodesids == ids(2))]; % Ordered nodes
    matPropind = truss.elems(i, end); % Material Prop
    x = truss.nodes(ids(1), :); % Node 1
    y = truss.nodes(ids(2), :); % Node 2
    xi = x(1); xj = y(1);
    length(i) = sqrt((xi - xj)^2);
    
    E = truss.mat(1, 1); % Youngs Modulus
    S = truss.mat(1, 2); % Surface Area
    
    % Local Stiffness
    ke =E*S*[1 -1;-1 1]/length(i); % Local K
    % Ordered Stiffness
    Kord(idsOrd,idsOrd) = Kord(idsOrd,idsOrd) + ke;
end

%% Loads

for n = 1:size(truss.loads, 1)
    ordFnode = find(ordNodesids==truss.loads(n,1));
    Ford(ordFnode) = truss.loads(n, 2); % X-Coordinate
end

if sol==1
    %% Method 1
    %---------------------
    % Boundary Conditions
    
    % Thined F and M
    bcremOrd = zeros(DoF, 1);
    for n = 1:size(truss.BC, 1)
      bcnode = truss.BC(n,1);
      ordBCnode = find(ordNodesids==bcnode);
      bcremOrd(ordBCnode) = truss.BC(n,2);
    end
    rmKord = Kord(~bcremOrd,~bcremOrd); % Removing the corresponding rows and columns of bcrem
    newFord = Ford(~bcremOrd); % Removing the corresponding rows of bcrem
    % New U as Matrix Solution to [K]{u} = {F}
    Uord = rmKord\newFord;
    
    % Rentering the previosuly removed nodal data, ordered
    j = 1;
    for i = 1:size(uord, 1)
      if bcremOrd(i) == 1
        uord(i) = 0; % 0 strain because of fixed support
      else
        uord(i) = Uord(j); % value from U as no fixed support
        j = j+1;
      end
    end
    for i=1:DoF
        uxyOrd(i,1) = uord(i);
    end
    
    % Decompose u
    ui = uxyOrd(1:truss.iinodes);
    ub = uxyOrd(truss.iinodes+1, end);
    
    %% Plot the Deformed Formed
    plottinOrd(uxyOrd); % Plot with Reordering
end
iinodes = truss.iinodes;
%% Decomposing Matrices
if last == 1
    iinodes = iinodes+1;
end

% Decompose f
fi = Ford(1:iinodes);
fb = Ford(iinodes+1:end);
fb = [-10e5;10e5];

% Decompose K
Kii = Kord(1:iinodes,1:iinodes);
Kib = Kord(1:iinodes,iinodes+1:end);
Kbi = Kord(iinodes+1:end, 1:iinodes);
Kbb = Kord(iinodes+1:end, iinodes+1:end);

%% Calculating maps
Sps = Kbb - (Kbi*(Kii\Kib));
bps = fb - (Kbi*(Kii\fi));
end