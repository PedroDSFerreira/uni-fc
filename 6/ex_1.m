clc
clear all
close all

% Dados
L = 50;
k = 0.93;
c = 0.094;
p = 8.9;

dx = 0.5;
dt = 0.1;
tf = 500;
Ti = 100;

% a) Método de Euler
t = 1:dt:tf;
x = 1:dx:L;
Nt = length(t);
Nx = length(x);

T = zeros(Nx,Nt);
T(2:end-1,1) = Ti;

eta = (k*dt)/(c*p*(dx^2));
for n=1:Nt-1
    for i=2:Nx-1
        T(i,n+1) = T(i,n) + eta*(T(i-1,n)-2*T(i,n)+T(i+1,n));
    end
end

index = find(x==L/4);
T_L4 = T(index,:);
T_50 = interp1(T_L4(index+2:end),t(index+2:end),50)
figure(3)
plot(t(index+2:end),T_L4(index+2:end))

figure(1)
mesh(t,x,T)
figure(2)
contour(x,t,T')