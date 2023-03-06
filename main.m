%% Add path
clear all; close all; clc
addpath('base/');
nbLocalElems = 10;
nbSub = 20;
plt = 1;
er_max = 1e-5;
max_iter = 100;


%% 1.1 FEM
nbSubFEM = 1;
[u, K, Kord, F, Ford, ~] = FEM(nbSubFEM, nbLocalElems, plt);
%% 2.1 Primal Schur - Direct
[u, u_BN, u_IN, ~, ~] = primalSchur(nbSub, nbLocalElems, plt);

%% 2.2 Conditioning of Direct Methods

% Primal
hH = [];
condSp = [];
hHK = [];
condK = [];
nbLocalElemsiter = [5:25];
nbLocalElemsiter = [50:70];

nbSubiter = [10:30];
for i=1:length(nbSubiter)
    [~, ~, ~, Sp, hi] = primalSchur(nbSubiter(i), nbLocalElemsiter(i), 0);
    [~, K, Kord, ~, ~, hk] = FEM(1, nbLocalElemsiter(i), 0);
    condSp = [condSp cond(Sp)];
    condK = [condK cond(K)];
    hH = [hH hi];
    hHK = [hHK hk];
end
figure(1)
plot(hH, condSp);
xlabel('h/H');
ylabel('Condition Number')
title('h/H Vs. Condition Number for Primal Schur');
saveas(figure(1), fullfile('assets/hHVSConditonPrim.png'));
figure(4)
plot(hHK, condK);
xlabel('h (m)');
ylabel('Condition Number')
title('h Vs. Condition Number for unstructured FEM');
saveas(figure(4), fullfile('assets/hVSConditonFEM.png'));

%% 2.3 Primal Schur - Conjugate Gradient
[u, iter, res] = conjGradPrimal(nbSub, nbLocalElems, 0);
figure(2)
plot([1:iter], res');
xlabel('Number of Iterations');
ylabel('Residual')
ylim([0.3,1.3]);
title('Residual Vs. Iterations for Conjugate Gradient Primal Schur');
saveas(figure(2), fullfile('assets/iterVSresidualPrimConjGrad.png'));

%% 2.4 Iterations taken for Conjugate Gradient
nbSubiter = [10:20];
iter = zeros(length(nbSubiter), 1);
for i=1:length(nbSubiter)
    [u, iter(i), ~] = conjGradPrimal(nbSubiter(i), 10, 0);
end

figure(3)
plot(nbSubiter, iter);
xlabel('Number of Subdomains');
ylabel('Iterations Taken')
title('Absence of Scalability in Conjugate Gradient Method for Primal Schur');
saveas(figure(3), fullfile('assets/iterVSnbSubConjPrimal.png'));

%% 2.5 Preconditioned Conjugate Gradient
[u, iter, Spinv, Sp] = conjGradPrePrim(nbSub, nbLocalElems, max_iter, er_max, plt);


%% 2.7 Dual Schur - Direct Method
[u, u_BN, u_IN, ~, ~] = dualSchur(nbSub, nbLocalElems, plt);

%% 2.8 FETI Method
[u, it, Spinv, Sp, Sd] = FETI(nbSub, nbLocalElems, max_iter, er_max, plt);

%% 2.9 LATIN Method
 [W_hat_concat, W_hat_concat_stack] = latin(nbSub, nbLocalElems, plt);