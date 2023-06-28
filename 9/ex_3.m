clc
clear all
close all

%   Condições iniciais
a = 1;
b = 4;
x_match = 0.3*a;
V0=10;              % Potencial poço

V = 0;              % [-a;a]
h = 1e-3;

x_prog = -b:h:x_match;
N_prog = length(x_prog);
y_prog = zeros(1,N_prog);
y_prog(1) = 0;
y_prog(2) = h;
V_prog = zeros(1,N_prog);

x_regr = x_match:h:b;
N_regr = length(x_regr);
y_regr = zeros(1,N_regr);
y_regr(end) = 0;
y_regr(end-1) = h;
V_regr = zeros(1,N_regr);


E_guess = 4*pi^2/8*[0.9,0.95];

tol = 1e-5;
max_iter = 1e4;
result = zeros(1,max_iter);

for i=1:N_prog
   if x_prog(i)<-a
      V_prog(i) = V0; 
   end
end
for i=1:N_regr
   if x_regr(i)>a
      V_regr(i) = V0; 
   end
end

for i=1:max_iter
    E=E_guess(i);
    
    g_prog=2*(E*ones(1,N_prog)-V_prog);
    
    for k=2:N_prog-1
       y_prog(k+1)= 1/(1+h^2/12*g_prog(k+1)) * (2*y_prog(k)*(1-5*h^2/12*g_prog(k)) -y_prog(k-1)*(1+h^2/12*g_prog(k-1)));
    end
    
    
    g_regr=2*(E*ones(1,N_regr)-V_regr);
    
    for k=N_regr-1:-1:2
       y_regr(k-1)= 1/(1+h^2/12*g_regr(k-1)) * (2*y_regr(k)*(1-5*h^2/12*g_regr(k)) -y_regr(k+1)*(1+h^2/12*g_regr(k+1)));
    end
    
    y_regr = y_regr./y_regr(1).*y_prog(end);
    
    dy_prog = 1/h*(25/12*y_prog(end)-4*y_prog(end-1)+3*y_prog(end-2)-4/3*y_prog(end-3)+1/4*y_prog(end-4));
    dy_regr = 1/h*(-25/12*y_regr(1)+4*y_regr(2)-3*y_regr(3)+4/3*y_regr(4)-1/4*y_regr(5));
    
    result(i) = ((dy_prog/y_prog(end))-(dy_regr/y_regr(1)))/((dy_prog/y_prog(end))+(dy_regr/y_regr(1)));
    
    if abs(result(i))<tol
        break
    end
    
    if i>1
       m =  (result(i)-result(i-1))/(E_guess(i)-E_guess(i-1));
       E_guess(i+1)=E_guess(i)- result(i)/m;
    end
end
E
plot(x_regr,y_regr)
hold on
plot(x_prog,y_prog)
