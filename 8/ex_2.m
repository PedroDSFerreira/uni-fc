clear all;
close all;
clc;

L = 1;
h = 0.025;

N = 2*L/h+1;

V_old = zeros(N,N);

n1 = (N-1)/4 + 1;
n2 = 3*(N-1)/4 + 1;

% inicializar parte oca de condensador
for i=n1:n2
    V_old(i,n1:n2) = 1;
end
V_new = V_old;

max_iter = 1e4;
tol = 1e-7;

for it = 1:max_iter
    for i=2:N-1
        for j=2:N-1
            if ~((i>=n1 || i<=n2) && j>=n1 && j<=n2) || ~((j>=n1 || j<=n2) && i>=n1 && i<=n2)
                V_new(i,j) = (V_old(i,j+1)+V_old(i,j-1)+V_old(i+1,j)+V_old(i-1,j))/4;    
            end
        end
    end
    % condição paragem
    if sqrt(sum(sum(V_new-V_old))).^2/sqrt(sum(sum(V_new))).^2 < tol
        n_iter = it;
        xy = -L:h:L;
        mesh(xy,xy,V_new)
        break
    end
    V_old=V_new;
end


[Ex,Ey]=gradient(V_new,h,h);
Ex=-Ex;
Ey=-Ey;
figure(3)
quiver(xy,xy,Ex,Ey)
axis equal


% Cálculo capacidade
epsilon_0 = 1;
E = sqrt(Ex(1,:).^2+Ey(1,:).^2); 
capacity = 4 * epsilon_0*trapz(E)*h
