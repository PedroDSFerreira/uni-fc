clear all
close all
clc

% Constantes
miu = 1e-3; %Kg/m
L = 1; %m
T = 1e3; %N
B = 0; %y(L)=0 -- é o objetivo (condição fronteira)
n = 1;

% Vetores
h = 0.01;
x = 0:h:L;
N = length(x);

A = full(gallery('tridiag',N,1,-2,1)); % criar matriz A

[eigen_vectors, eigen_values]=eigs(A,3,'sm');

for i = 1:length(eigen_values)
    w(i) = sqrt(-eigen_values(i,i)*T/(miu*h^2));
end

%%só para confirmar se dá inteiros
for i = 1:length(eigen_values)
    nn(i) = w(i)/w(1);
end

plot(x,eigen_vectors(:,2))