function [x,v] = Runge_Kutta_2(x, v, t, N, h, col, matriz, pesos, arg_equacoes)
    
    % equações necessárias
    dxdt = @(T,X,V) arg_equacoes;
    dvdt = @(T,X,V) arg_equacoes;

    for i=1:N-1

        % primeira derivação
        r1x = dxdt(t(i),x(i),v(i));     
        r1v = dvdt(t(i),x(i),v(i));

        % estimativa de valores no centro do intervalo
        x2 = x(i) + r1x*matriz(2,1)*h;
        v2 = v(i) + r1v*matriz(2,1)*h;
        t2 = t(i) + col(2)*h;
        
        % segunda derivação
        r2x = dxdt(t2,x2,v2);
        r2v = dvdt(t2,x2,v2);


        x(i+1) = x(i) + (pesos(1)*r1x + pesos(2)*r2x)*h;
        v(i+1) = v(i) + (pesos(1)*r1v + pesos(2)*r2v)*h;
    end

end