% Analytical solution 
clear all; close all;
% Variables
addpath('base/');
run('data.m'); % Generates truss mesh
Fd = 10e5; % Force on the end node

E = truss.mat(1, 1); % Youngs Modulus
S = truss.mat(1, 2); % Surface Area

x = linspace(0,100,10);
dx = x.*(Fd/(S*E));
figure
plot(x,dx)
title("Tensile Response of a Beam")
xlabel("x")
ylabel("dx")
