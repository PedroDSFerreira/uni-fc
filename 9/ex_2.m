clc
close all
clear all

%   Condições iniciais
h=1e-2;         %passo
r_max=50;
r=0:h:r_max;    %raios do átomo de hidrogénio
N=length(r);

u=zeros(1,N);
u(end)=0;       %condição fronteira
u(end-1)=h;

n=3;
l=1;            % valores possíveis: l=[0,n-1]

% E = -1/2*n^2;
result=zeros(1,N);

tol=1e-5;       %tolerância de valores
max_iter=1e4;   %iteração máxima

E_guess=-1/20*[0.9,0.95];

for i=1:max_iter
    
    E=E_guess(i);
	g=2*(E-(l*(l+1))./(2*r.^2)+1./r);
    
    for k=N-1:-1:3
       u(k-1)= 1/(1+h^2/12*g(k-1)) * (2*u(k)*(1-5*h^2/12*g(k)) -u(k+1)*(1+h^2/12*g(k+1)));
    end
    u(1) = interp1(r(2:5),u(2:5),0,'spline');
    
    C=trapz(u.^2)*h;
    u=u/sqrt(C);  % normalização
    
    % shooting
    result(i) = u(1);
    
    if abs(result(i))<tol
        break
    end
    
    if i>1
       m =  (result(i)-result(i-1))/(E_guess(i)-E_guess(i-1));
       E_guess(i+1)=E_guess(i)- result(i)/m;
    end

end
R=u./r;
plot(r,R)
E
figure()
plot(r,u)