clc
close all
clear all

%   Condições iniciais
a=1;
V=0;     % ]-a;a[
h=1e-2;
x=-a:h:a;
N=length(x);
y=zeros(1,N);
y(2)=h;

E_guess=pi^2/8*[0.9,1.1];
result=zeros(1,N);
tol=1e-5;
max_iter=1e4;

for i=1:max_iter
    
    E=E_guess(i);
	g=2*E*ones(1,N);
    
    for k=2:N-1
       y(k+1)= 1/(1+h^2/12*g(k+1)) * (2*y(k)*(1-5*h^2/12*g(k)) -y(k-1)*(1+h^2/12*g(k-1)));
    end
    
    C=trapz(y.^2)*h;
    y=y/sqrt(C);  % normalização
    
    % shooting
    result(i) = y(end);
    
    if abs(result(i))<tol
        break
    end
    
    if i>1
       m =  result(i)-result(i-1)/(E_guess(i)-E_guess(i-1));
       E_guess(i+1)=E_guess(i)- result(i)/m;
    end

end

E
plot(x,y)