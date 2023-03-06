function [W_hat_concat, W_hat_concat_stack] = latin(nbSub, nbLocalElem, plt)
    % Solver for monoscale latin approach
    %  [W_hat_concat, W_hat_concat_stack] = latin(nbSub, nbLocalElem)
    % Output :  W_hat_concat = final displacement vector
    %           W_hat_concat_stack = displacement vector at each iteration
    % Input :   nbSub = number of subdomains
    %           nbLocalElems = number of elements in subdomain
    %           plt = 1 to plot figures, 0 otherwise
    % Example : [W_hat_concat, W_hat_concat_stack] = latin(10, 1)

    %% Initialization of the problem
    truss = mesher(nbSub, nbLocalElem);
    iter = truss.nbSub;
    Fd = truss.Fd; % Force on the end node


    %DoF = truss.Dim*truss.nbNodes;
    DoF = 2; % 2 for 1D elements - only two boundary terms per element
    E = truss.mat(1, 1); % Youngs Modulus
    S = truss.mat(1, 2); % Surface Area
    L = truss.L; % Length of the element
    ko = E*S/L;
    ke = ko*[1 -1;-1 1]; % Local K


    % Interface Stiffnesses
    nbInterfaces = truss.nbSub + 1;
    InterfaceStiff = zeros(nbInterfaces,1);
    InterfaceStiff(1) = 1000*ko; % Infinite Stiffness - rigid surface
    for i = 2:nbInterfaces-1
        InterfaceStiff(i) = ko/(i-1);
    end

    %Concatenated Interface Stiffness matrices
    Kinv_concat = [];
    for i = 1:truss.nbSub
        Kint = zeros(size(ke)); 
        Kint(1,1) = InterfaceStiff(i);
        Kint(end,end) = InterfaceStiff(i+1);
        Kinv_concat = [Kinv_concat inv(Kint+ke)];
    end

    %% Initialization 
    % Displacement
    W_concat = zeros(DoF*truss.nbSub,1);
    W_hat_concat = zeros(nbInterfaces,1);
    W_hat_concat_stack = W_hat_concat';
    % force
    F_concat = zeros(DoF*truss.nbSub,1);
    F_hat_concat = zeros(nbInterfaces,1);
    F_hat_concat_stack = F_hat_concat';
    %% Loop on Number of Iterations
    for it = 1:iter
        % Local Stage, First Iteration 
        % Border - Imposed Force
        F_hat_concat(end) = Fd;
        % Interior Elements - Imposed Displacement 
        for i = 1:truss.nbSub
            W_hat_concat(i) = W_concat(DoF*i-1); 
        end
        
        % Border - Calculated Displacement
        W_hat_concat(end) = W_concat(end) + InterfaceStiff(1)*(F_hat_concat(end)-F_concat(end));
        % Interior Elements - Calculated Force
        for i = 1:truss.nbSub
            F_hat_concat(i) = -F_concat(DoF*i-1) + InterfaceStiff(i)*(W_hat_concat(i)-W_concat(DoF*i-1)); 
        end
        LatinRHS = zeros(DoF,1);
        % Linear Stage, first Iteration 
        for i = 1:truss.nbSub
            LatinRHS(1) = -F_hat_concat(i) + InterfaceStiff(i)*W_hat_concat(i);
            LatinRHS(end) = F_hat_concat(i+1) + InterfaceStiff(i+1)*W_hat_concat(i+1);
            W_concat(DoF*i-1:DoF*i) = Kinv_concat(:,DoF*i-1:DoF*i)*LatinRHS;
            F_concat(DoF*i) = F_hat_concat(i+1) + InterfaceStiff(i+1)*(W_hat_concat(i+1)-W_concat(DoF*i)); 
            F_concat(DoF*i-1) = -F_hat_concat(i) + InterfaceStiff(i)*(W_hat_concat(i)-W_concat(DoF*i-1)); 
        end
        W_hat_concat_stack = [W_hat_concat_stack;W_hat_concat'];
        F_hat_concat_stack = [F_hat_concat_stack;F_hat_concat'];
    end

    % Plotting
    if plt == 1
        % coordinates
        x_coor = linspace(0,truss.nbSub*L,truss.nbSub+1)';
        figure(1);plot(x_coor,W_hat_concat);xlabel("X (m)");ylabel("Displacement (m)");title("Latin Method Displacements");
        figure(2);spy(W_hat_concat_stack);xlabel("boundary terms");ylabel("iterations");title("Nonzero Boundary Displacement terms");
        figure(3);spy(F_hat_concat_stack);xlabel("boundary terms");ylabel("iterations");title("Nonzero Boundary Force terms");

        % save figures
        saveas(figure(1), fullfile('assets/latin_displacements.png'));
        saveas(figure(2), fullfile('assets/latin_displacements_spy.png'));
        saveas(figure(3), fullfile('assets/latin_forces_spy.png'));
    end
end
