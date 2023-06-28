%% Executar antes de cada secção

clc;
close all;
clear all;

b = [0 ; 0; 0 ; 0; 0 ; 1000; 0 ; 0 ];

tol = 1e-7;

A = [  -1   0    0   sqrt(2)/2   1    0     0        0;
        0   -1   0   sqrt(2)/2   0    0     0        0;
        0   0   -1       0       0    0    0.5       0;
        0   0    0  -sqrt(2)/2   0   -1    0.5       0;
        0   0    0       0      -1    0     0        1;
        0   0    0       0       0    1     0        0;
        0   0    0  -sqrt(2)/2   0    0  sqrt(3)/2   0;
        0   0    0       0       0    0 -sqrt(3)/2  -1];
iter_max = 1000;
size = 8;
x_old = ones(size,1);
x_new = x_old;

%% Método de Jacobi

for it = 1:iter_max
    for i=1:size
        x_new(i)= b(i)/A(i,i);    
        for j = 1:size
            if i~= j
                x_new(i)= x_new(i)-1/A(i,i)*A(i,j)*x_old(j);
            end
        end
    end
    if max(abs(x_new-x_old))/max(x_new) < tol
        it
        x_new
        break
    end
    x_old = x_new;
end
'Método não convergiu'
%% Método de Gauss-Seidel

for it = 1:iter_max
    for i=1:size
        x_new(i)= b(i)/A(i,i);    
        for j = 1:size
            if i~= j
                x_new(i)= x_new(i)-1/A(i,i)*A(i,j)*x_new(j);
            end
        end
    end
    if max(abs(x_new-x_old))/max(x_new) < tol
        it
        x_new
        break
    end
    x_old = x_new;
end
'Método não convergiu'
%% Método de sub- ou sobre-relaxação

alpha = 1.25;

for it = 1:iter_max
    for i=1:size
        x_new(i)= (1-alpha)*x_old(i);    
        for j = 1:size
            if i~= j
                x_new(i)= x_new(i)-alpha/A(i,i)*(A(i,j)*x_new(j)+b(i));
            end
        end
    end
    if max(abs(x_new-x_old))/max(x_new) < tol
        it
        x_new
        break
    end
    x_old = x_new;
end
'Método não convergiu'