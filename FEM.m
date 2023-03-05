function [u, K, Kord, F, Ford, h] = FEM(nbSub, nbLocalElems, plt)

truss = mesher(nbSub, nbLocalElems);
h = truss.h;% Generate the Truss
%% Initialize Matrices

DoF = truss.Dim*truss.nbNodes;
K = sparse(DoF, DoF);
Kord = sparse(DoF, DoF);
length = zeros(truss.nbElems,1);
u = sparse(DoF,1);
uord = sparse(DoF,1);
F = sparse(DoF,1);
Ford = sparse(DoF,1);

%% Local Calculations

for i=1:truss.nbElems
    ids = truss.elems(i, (1:end-1)); % Node IDs
    idsOrd = [find(truss.ordNodesids == ids(1)) find(truss.ordNodesids == ids(2))]; % Ordered nodes
    matPropind = truss.elems(i, end); % Material Prop
    x = truss.nodes(ids(1), :); % Node 1
    y = truss.nodes(ids(2), :); % Node 2
    xi = x(1); xj = y(1);
    length(i) = sqrt((xi - xj)^2);
    
    E = truss.mat(1, 1); % Youngs Modulus
    S = truss.mat(1, 2); % Surface Area
    
    % Local Stiffness
    ke =E*S*[1 -1;-1 1]/length(i); % Local K
    
    % Original Stiffness
    K(ids,ids) = K(ids,ids) + ke;

    % Ordered Stiffness
    Kord(idsOrd,idsOrd) = Kord(idsOrd,idsOrd) + ke;
end


%% Loads

for n = 1:size(truss.loads, 1)
    F(truss.loads(n,1)) = truss.loads(n, 2); % X-Coordinate
    ordFnode = find(truss.ordNodesids==truss.loads(n,1));
    Ford(ordFnode) = truss.loads(n, 2); % X-Coordinate
end

%% Method 1
%---------------------
% Boundary Conditions

% Thined F and M
bcrem = zeros(DoF, 1);
bcremOrd = zeros(DoF, 1);
for n = 1:size(truss.BC, 1)
  bcnode = truss.BC(n,1);
  ordBCnode = find(truss.ordNodesids==bcnode);
  bcrem(bcnode) = truss.BC(n,2);
  bcremOrd(ordBCnode) = truss.BC(n,2);
end

rmK = K(~bcrem, ~bcrem); % Removing the corresponding rows and columns of bcrem
newF = F(~bcrem); % Removing the corresponding rows of bcrem
rmKord = Kord(~bcremOrd,~bcremOrd);
newFord = Ford(~bcremOrd);
% New U as Matrix Solution to [K]{u} = {F}
U = rmK\newF;
Uord = rmKord\newFord;

% Rentering the previosuly removed nodal data
j = 1;
for i = 1:size(u, 1)
  if bcrem(i) == 1
    u(i) = 0; % 0 strain because of fixed support
  else
    u(i) = U(j); % value from U as no fixed support
    j = j+1;
  end
end
for i=1:DoF
    uxy(i,1) = u(i);
end

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
if plt == 1
figure
spy(K);
title('Sparsity Patterns for K');
saveas(figure(1), fullfile('assets/sparse_K.png'));
figure
spy(Kord);
title('Sparsity Patterns for K_{ord}');
saveas(figure(2), fullfile('assets/sparse_Kord.png'));

%% Plot the Deformed Formed
if nbSub == 1;
plottin(truss,uxy)
title('Displacement without Reordering')% Plot without Reordering
saveas(figure(3), fullfile('assets/fem_uxy.png'));

plottinOrd(truss, uxyOrd) % Plot with Reordering
saveas(figure(4), fullfile('assets/fem_uxy.png'));

end
end
end