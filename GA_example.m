clc
clear all
tic
lb = [zeros(1,10),0.0001];
ub = [1000 1000 1 1000 1000 1000 1000 1000 1000 10 1];
% x =[landa   nu gamma vu   g0    g1    g2    g3 epstar S0 kappa]
% x=[100 500 0.5 0.01 700 700 600 10 0.5 0.5 500]

% landa=x(1);
% nu=x(2);
% gamma=x(3);
% vu=x(4);
% g0=x(5);
% g1=x(6);
% g2=x(7);
% g3=x(8);
% epstar=x(9);
% S0=x(10);
% kappa=x(11)

FitnessFunction = @(x) GA_init(x);
options = optimoptions('ga',...
    'PlotFcns',  {@gaplotbestf @gaplotbestindiv },'MaxGenerations',10,...
    'FunctionTolerance',1e-6);

numberOfVariables = 11;
[x,exitflag,output] = ga(FitnessFunction,numberOfVariables,[],[],[],[],lb,ub,[],options)
toc
