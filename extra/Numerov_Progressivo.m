function [y] = Numerov_Progressivo(y, g, S, N, h)
    for k=2:N-1
       y(k+1)= 1/(1+h^2/12*g(k+1)) * (2*y(k)*(1-5*h^2/12*g(k)) -y(k-1)*(1+h^2/12*g(k-1)) +h^2/12*(S(k-1)+10*S(k)+S(k+1)));
    end
end
