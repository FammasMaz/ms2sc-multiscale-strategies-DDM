function [x, dx] = solExact(nbSub, nbLocalElems)
    % Solves the exact solution for a truss with a force applied to the end node
    % Inputs:   nbSub: Number of subdomains
    %           nbLocalElems: Number of elements per subdomain
    % Outputs:  x: Vector of position of the nodes
    %           dx: Vector of displacement of the nodes



    truss = mesher(nbSub, nbLocalElems); % Generates truss mesh
    Fd = truss.Fd; % Force on the end node

    E = truss.mat(1, 1); % Youngs Modulus
    S = truss.mat(1, 2); % Surface Area

    x = linspace(0,100,10);
    dx = x.*(Fd/(S*E));
    figure
    plot(x,dx)
    title("Tensile Response of a Beam")
    xlabel("X (m)")
    ylabel("Displacement (m)")
    saveas(figure(1), fullfile("assets/Tensile Response of a Beam.png"));
end