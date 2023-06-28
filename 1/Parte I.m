clc
clear all
close all

%% Queda de um grave
% Constantes

m = 0.15; %kg
g = 9.8;  %kg/s^2

%
dt = 0.2;
ti = 0;
tf = 10;

t = ti:dt:tf;
v = zeros(1,length(t));
z = zeros(1,length(t));

v(1) = 0;
z(1) = 6;

count = 1;

while (count ~= length(t))
   v(count+1) = v(count) -g*dt;
   z(count+1) = z(count) +v(count)*dt; 
   count = count+1;
end

plot(t,v)
hold on
plot(t,v(1)-g*t, 'x')

figure()
hold on
plot(t,z)
plot(t,z(1)-(g/2)*t.^2)
t_= interp1(z(2:end), t(2:end),0);

%% Oscilador harmónico simples
%Constantes
K = 2.5; %N/m
m = 0.5; %kg

%
dt = 0.1;
ti = 0;
tf = 10;

t = ti:dt:tf;
v = zeros(1,length(t));
z = zeros(1,length(t));

v(1) = 0;
z(1) = 0.1;

count = 1;

while (count ~= length(t))
   v(count+1) = v(count) -(K/m)*z(count)*dt;
   z(count+1) = z(count) +v(count)*dt; 
   count = count+1;
end

figure()
plot(t,z)
hold on
plot(t,z(1)*cos(sqrt(K/m)*t))
figure()
plot(t,v)
hold on
plot(t,-z(1)*sqrt(K/m)*sin(sqrt(K/m)*t))
