function [x,v] = Runge_Kutta_4(x, v, t, N, h, col, matriz, pesos, arg_equacoes)

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
        
        % estimativa de valores no centro do intervalo
        x3 = x(i) + (r1x*matriz(3,1)+r2x*matriz(3,2))*h;
        v3 = v(i) + (r1v*matriz(3,1)+r2v*matriz(3,2))*h;
        t3 = t(i) + col(3)*h;

        % terceira derivação
        r3x = dxdt(t3,x3,v3);
        r3v = dvdt(t3,x3,v3);

        % estimativa de valores no centro do intervalo
        x4 = x(i) + (r1x*matriz(4,1)+r2x*matriz(4,2)+r3x*matriz(4,3))*h;
        v4 = v(i) + (r1v*matriz(4,1)+r2v*matriz(4,2)+r3v*matriz(4,3))*h;
        t4 = t(i) + col(4)*h;
        
        % quarta derivação
        r4x = dxdt(t4,x4,v4);
        r4v = dvdt(t4,x4,v4);


        x(i+1) = x(i) + (pesos(1)*r1x + pesos(2)*r2x + pesos(3)*r3x + pesos(4)*r4x)*h;
        v(i+1) = v(i) + (pesos(1)*r1v + pesos(2)*r2v + pesos(3)*r3v + pesos(4)*r4v)*h;
    end

end