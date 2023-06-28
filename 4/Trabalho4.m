clear all
close all
clc

%% 4.1 - Runge Kutta (3º)

% Dados
EPSILON = 0.1;
tf = 100;
h = 0.01;
v0 = 0.7;
x0 = 0.2;

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;

[x, v, t] = runge_kutta3(x,v,t,h,N,EPSILON);

% a)
figure(1)

plot(t,x)
title("4.1 a) (y/t)")
xlabel('t')
ylabel('y')


figure(2)

plot(x,v)
title("4.1 a) (v/y)")
xlabel('y')
ylabel('v')

% b)

% Dados
EPSILON = 1;
v0 = 7;
x0 = 2;

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;

[x, v, t] = runge_kutta3(x,v,t,h,N,EPSILON);

figure(3)

plot(t,x)
title("4.1 b) (y/t)")
xlabel('t')
ylabel('y')


figure(4)

plot(x,v)
title("4.1 b) (v/y)")
xlabel('y')
ylabel('v')


%% 4.2 - Runge Kutta (4º)

% a)
% Dados
EPSILON = 0.3;
tf = 100;
h = 0.01;
v0 = 2;
x0 = 0;
F0 = 1;

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;

[x, v, t] = runge_kutta4(x,v,t,h,N,EPSILON,F0);

figure(4)

plot(t,x)
title("4.2 a) (y/t)")
xlabel('t')
ylabel('y')


figure(5)

plot(x,v)
title("4.2 a) (v/y)")
xlabel('y')
ylabel('v')

% b)
F0 = 1.5;

t = 0:h:tf;
N = length(t);
x = zeros(1,N);
v = zeros(1,N);

x(1) = x0;
v(1) = v0;

[x, v, t] = runge_kutta4(x,v,t,h,N,EPSILON,F0);

figure(6)

plot(t,x)
title("4.2 b) (y/t)")
xlabel('t')
ylabel('y')


figure(7)

plot(x,v)
title("4.2 b) (v/y)")
xlabel('y')
ylabel('v')


%% 4.3 - ode45
EPSILON = 0.1;
f1 = @(T,XV) [XV(2); -XV(1)-EPSILON*(XV(1)^2-1)*XV(2)];

EPSILON = 1;
f2 = @(T,XV) [XV(2); -XV(1)-EPSILON*(XV(1)^2-1)*XV(2)];


tspan = [0 100];

reltol = 3E-14;
abstol_1 = 1E-13;
abstol_2 = 1E-12;
options = odeset('RelTol',reltol,'AbsTol',[abstol_1 abstol_2]);

%condições iniciais
y0_1 = 1;
y0_2 = 8;
y0_3 = -5;

v0_1 = 7;
v0_2 = 4;
v0_3 = 2;


% epsilon=0.1
[t1,xv1] = ode45(f1, tspan, [y0_1 v0_1], options);
[t2,xv2] = ode45(f1, tspan, [y0_2 v0_2], options);
[t3,xv3] = ode45(f1, tspan, [y0_3 v0_3], options);


% epsilon=1
[t4,xv4] = ode45(f2, tspan, [y0_1 v0_1], options);
[t5,xv5] = ode45(f2, tspan, [y0_2 v0_2], options);
[t6,xv6] = ode45(f2, tspan, [y0_3 v0_3], options); 


figure(8)
hold on
plot(xv1(:,1),xv1(:,2))
plot(xv2(:,1),xv2(:,2))
plot(xv3(:,1),xv3(:,2))

title("4.3  e=0.1")
xlabel('y')
ylabel('v')


figure(9)
hold on
plot(xv4(:,1),xv4(:,2))
plot(xv5(:,1),xv5(:,2))
plot(xv6(:,1),xv6(:,2))

title("4.3  e=1")
xlabel('y')
ylabel('v')


%% Funções utilizadas

function [x,v,t] = runge_kutta3(x, v, t, h, N, EPSILON)
    dxdt = @(T,X,V) V;
    dvdt = @(T,X,V) -X-EPSILON*(X^2-1)*V ;

    for i=1:N-1

        % primeira derivação
        r1x = dxdt(t(i),x(i),v(i));
        r1v = dvdt(t(i),x(i),v(i));

        % estimativa de valores no centro do intervalo
        x2 = x(i) + r1x*h;
        v2 = v(i) + r1v*h;
        t2 = t(i) + h;

        % segunda derivação
        r2x = dxdt(t2,x2,v2);
        r2v = dvdt(t2,x2,v2);

        x3 = x(i) + r1x*h/4 + r2x*h/4;
        v3 = v(i) + r1v*h/4 + r2v*h/4;
        t3 = t(i) + h/2;

        % terceira derivação
        r3x = dxdt(t3,x3,v3);
        r3v = dvdt(t3,x3,v3);


        x(i+1) = x(i) + (1*r1x + 1*r2x + 4*r3x)*h/6;
        v(i+1) = v(i) + (1*r1v + 1*r2v + 4*r3v)*h/6;

    end
end

function [x,v,t] = runge_kutta4(x, v, t, h, N, EPSILON, F0)

    dxdt = @(T,X,V) V;
    dvdt = @(T,X,V) -X-EPSILON*(X^2-1)*V+F0*cos(1.7*T) ;

    for i=1:N-1

        % primeira derivação
        r1x = dxdt(t(i),x(i),v(i));
        r1v = dvdt(t(i),x(i),v(i));

        x2 = x(i) + r1x*h/3;
        v2 = v(i) + r1v*h/3;
        t2 = t(i) + h/3;

        % segunda derivação
        r2x = dxdt(t2,x2,v2);
        r2v = dvdt(t2,x2,v2);

        x3 = x(i) - r1x*h/3 + r2x*h;
        v3 = v(i) - r1v*h/3 + r2v*h;
        t3 = t(i) + h*2/3;

        % terceira derivação
        r3x = dxdt(t3,x3,v3);
        r3v = dvdt(t3,x3,v3);
        
        x4 = x(i) + r1x* - r2x*h + r3x*h;
        v4 = v(i) + r1v* - r2v*h + r3v*h;
        t4 = t(i) + h;

        % quarta derivação
        r4x = dxdt(t4,x4,v4);
        r4v = dvdt(t4,x4,v4);


        x(i+1) = x(i) + (1*r1x + 3*r2x + 3*r3x + 1*r4x)*h/8;
        v(i+1) = v(i) + (1*r1v + 3*r2v + 3*r3v + 1*r4v)*h/8;

    end
end