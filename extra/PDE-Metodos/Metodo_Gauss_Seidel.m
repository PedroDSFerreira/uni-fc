function [x] = Metodo_Gauss_Seidel(A, b, N, iter_max, tol)
    x_old = ones(N,1);
    x_new = x_old;
    
    for i = 1:iter_max
        for j=1:N
            x_new(j)= b(j)/A(j,j);    
            for k = 1:N
                if j~= k
                    x_new(j)= x_new(j)-1/A(j,j)*A(j,k)*x_new(k);
                end
            end
        end
        if max(abs(x_new-x_old))/max(x_new) < tol
            x = x_new;
            break
        end
        x_old = x_new;
    end
end