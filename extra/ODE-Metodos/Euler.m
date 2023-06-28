function [x,v] = Euler(x, v, t, N, arg_equacoes)

    % equações necessárias
    dvdt = @(T,X,V) arg_equacoes;

    for i=1:N-1
        v(i+1) = v(i) + dvdt(t(i),x(i),v(i))*h;
        x(i+1) = x(i) + v(i)*h;
    end
end