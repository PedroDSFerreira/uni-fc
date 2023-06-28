clc
clear all
close all

%% Dados

h = 0.1;
tf = 10;

K = 16;
m = 1;
x0 = 1;
v0 = 0;

%% Runge Kutta (2º)

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;

for i=1:N-1
    
    % primeira derivação
    r1x = v(i);     
    r1v = -(K/m)*x(i);
    
    % estimativa de valores no centro do intervalo
    x2 = x(i) + r1x*h/2;
    v2 = v(i) + r1v*h/2;
    
    % segunda derivação
    r2x = v2;
    r2v = -(K/m)*x2;
    

    x(i+1) = x(i) + (0*r1x + 1*r2x)*h;
    v(i+1) = v(i) + (0*r1v + 1*r2v)*h;
    
    
    
%     plot(t(i), x(i), 'o')
%     grid on
%     axis([-1 5 -2 2])
%     hold off
%     pause(0.5)
    
end


%% Runge Kutta (2º, simplificação)                  <----------------------
dxdt = @(T,X,V) V;
dvdt = @(T,X,V) -(K/m)*X;

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;


for i=1:N-1
    
    % primeira derivação
    r1x = dxdt(t(i),x(i),v(i));     
    r1v = dvdt(t(i),x(i),v(i));
    
    % estimativa de valores no centro do intervalo
    x2 = x(i) + r1x*h/2;
    v2 = v(i) + r1v*h/2;
    
    t2 = t(i) + h/2;
    
    % segunda derivação
    r2x = dxdt(t2,x2,v2);
    r2v = dvdt(t2,x2,v2);
    

    x(i+1) = x(i) + (0*r1x + 1*r2x)*h;
    v(i+1) = v(i) + (0*r1v + 1*r2v)*h;
end

%% Runge Kutta (2º, simplificação)                  <----------------------
dxdt = @(V) V;
dvdt = @(X) -(K/m)*X;

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;


for i=1:N-1
    
    % primeira derivação
    r1x = dxdt(v(i));     
    r1v = dvdt(x(i));

    % segunda derivação
    r2x = dxdt(v(i) + r1v*h/2);
    r2v = dvdt(x(i) + r1x*h/2);
    

    x(i+1) = x(i) + (0*r1x + 1*r2x)*h;
    v(i+1) = v(i) + (0*r1v + 1*r2v)*h;
end

%% Runge Kutta (4º, simplificação)
dxdt = @(T,X,V) V;
dvdt = @(X) -(K/m)*X;

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;


for i=1:N-1
    
    % primeira derivação
    r1x = dxdt(v(i));     
    r1v = dvdt(x(i));

    % segunda derivação
    r2x = dxdt(v(i) + r1v*h/2);
    r2v = dvdt(x(i) + r1x*h/2);
    
    % terceira derivação
    r3x = dxdt(v(i) + r2v*h/2);
    r3v = dvdt(x(i) + r2x*h/2);
    
    % quarta derivação
    r4x = dxdt(v(i) + r3v*h/2);
    r4v = dvdt(x(i) + r3x*h/2);
    
    %coeficientes errados!!!
    x(i+1) = x(i) + (0*r1x + 1*r2x + 2*r3x + 3*r4x)*h;
    v(i+1) = v(i) + (0*r1v + 1*r2v + 2*r3v + 3*r4v)*h;
end

