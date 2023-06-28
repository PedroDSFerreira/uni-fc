function [x,v] = Crank_Nicolson(x, v, t, N, A, arg_equacoes)
    
    % v(k+1)=v(k)+[f(t(k),y(k))+f(t(k+1),y(k+1))]*(h/2)
    % y(k+1)=y(k)+(v(k)+v(k+1))*(h/2)

    for i=1:N-1
        b = [1;1]; %mudar de acordo com problema,  -->AZ=b
        Z = linsolve(A,b);
        x(i+1) = Z(1);
        v(i+1) = Z(2);
    end
end