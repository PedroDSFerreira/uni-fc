close all
clear all
clc

% Dados
L = 50;
k = 0.93;
c = 0.094;
p = 8.9;

dx = 0.5;
dt = 0.1;
tf = 500;
Ti = 100;

% Método de Crank-Nicolson
t = 1:dt:tf;
x = 1:dx:L;

Nt = length(t);
Nx = length(x);

T = zeros(Nx,Nt);
T(2:end-1,1) = Ti;

eta = (k*dt)/(c*p*(dx^2));

b = zeros(Nx,1);
A = full(gallery('tridiag',Nx,-1,(2/eta)+2,-1)); % criar matriz tridiagonal


for n=1:Nt-1
    for i=2:Nx-1
        b(i) = T(i,n)+(2/eta-2)*T(i+1,n)+T(i+2,n);
    end
    b(1) = b(1)+T(1,n+1);
    b(Nx-2) = b(Nx-2)+T(Nx,n+1);
    Z = linsolve(A,b);
    T(2:Nx-1,n+1)=Z(2:Nx-1);
end

mesh(t,x,T)