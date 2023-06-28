clc
clear all
close all

%% Movimento a duas dimensões
% Constantes

m = 1; %kg
g = 9.8;  %kg/s^2
v0 = 50;
h0 = 0;
ALPHA = 37;

%

dt = 0.01;
ti = 0;
tf = 10;


t = ti:dt:tf;
v = zeros(2,length(t));
x_y = zeros(2,length(t));

v(:,1) = [v0*cos(ALPHA*pi/180);v0*sin(ALPHA*pi/180)];
x_y(:,1) = [0;h0];

count = 1;

while (x_y(2,count) >= 0 )
   v(1,count+1) = v(1,count);
   v(2,count+1) = v(2,count) - g*dt;
   x_y(1,count+1) = x_y(1,count) + v(1,count)*dt;
   x_y(2,count+1) = x_y(2,count) + v(2,count)*dt; 

   count = count+1;
end

x_y = x_y(:,2:count);
t = t(2:count);
v = v(2:count);

plot(x_y(1,:), x_y(2,:), 'o')

tvoo = interp1(x_y(2,:), t, 0)
alcance = interp1(x_y(2,:), x_y(1,:), 0)

R = abs(((v0^2)/g)*sin(2*ALPHA))


%% Movimento a Três Dimensões

% Constantes
ALPHA=1;
dt = 0.01;
ti = 0;
tf = 10;
g0 = [0; 0; -9.8];
OMEGA = [-7.292*10^(-5)*cos(ALPHA); 0; 7.292*10^(-5)*sin(ALPHA)];
R = [0; 0; 6.37*10^6];
H = [0; 0; 200];
v0 = [0; 0; 0];
m = 1;

t = ti:dt:tf;

axis = zeros(3,length(t));
v = zeros(3,length(t));

axis(:,1) = [0; 0; H(3)];
v(:,1) = [0; 0; 0];

for i=1:length(t)-1
    v(:,i+1) = v(:,i) + (g0 - 2*cross(OMEGA,v(:,i)) - cross(OMEGA,cross(OMEGA,R)))*dt;
    axis(:,i+1) = axis(:,i) + v(:,i)*dt;
    
    if axis(3,i+1)<0
        axis = axis(:,1:i+1);
        v = v(:,1:i+1);
        t = t(1:i+1);
        break
    end
end
axis(3,:)
plot3(axis(1,:), axis(2,:), axis(3,:))
t_emb = interp1(axis(3,end-1:end),t(end-1:end),0)
