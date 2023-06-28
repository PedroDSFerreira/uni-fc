close all
clear all
clc

%%

miu =1e-3;
L = 1;
T = 1e3;
n = 1;
omega = n*pi/L*sqrt(T/miu);

h = 0.01;

x = 0:h:L;
N = length(x);
y = zeros(1,N);
dydx = zeros(1,N);

y(1) = 0;
dydx(1) = 1;


d2ydt2 = @(Y) -omega^2*miu*Y/T;


% metodo euler-cromer
for i=1:N-1
    dydx(i+1) = dydx(i) + d2ydt2(y(i))*h;
    y(i+1) = y(i) + dydx(i+1)*h;
end
omega
figure()
plot(x,y,'b')

hold on

%% Método shooting
clear all

miu =1e-3;
L = 1;
T = 1e3;
n = 1;
h = 0.01;
x = 0:h:L;
N = length(x);
y = zeros(1,N);
dydx = zeros(1,N);

y(1) = 0;
dydx(1) = 1;

result = zeros(1,N);
guess = zeros(1,N);
guess(1) = 2000;
tol = 1e-7;     %tolerância

for j=1:N-1
    w = guess(j);
    d2ydt2 = @(Y) -w^2*miu*Y/T;
    
    for i=1:N-1
        dydx(i+1) = dydx(i) + d2ydt2(y(i))*h;
        y(i+1) = y(i) + dydx(i+1)*h;
    end
    result(j) = y(end);
    
    if abs(result(j))<tol
        guess = guess(1:j);
        result = result(1:j);

        break
    end
    
    plot(x,y,'-.g')
    
    if j>1
        m = (result(j)-result(j-1))/(guess(j)-guess(j-1));
        b = result(j)-m*guess(j);
        guess(j+1) = guess(j)-result(j)/m;
    end
end

omega_result = guess(end)
plot(x,y,'xr')

