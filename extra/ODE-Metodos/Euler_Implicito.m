function [x,v] = Euler_Implicito(x, v, t, N, A, arg_equacoes)
    
    % 𝑣𝑘+1 = 𝑣𝑘 + 𝑓(𝑡𝑘+1, 𝑦𝑘+1, 𝑣𝑘+1 ) ∗ ℎ
    % 𝑦𝑘+1 = 𝑦𝑘 + 𝑣𝑘+1 ∗ ℎ

    for i=1:N-1
        b = [arg_equacoes;arg_equacoes]; %mudar de acordo com problema,  -->AZ=b
        Z = linsolve(A,b);
        x(i+1) = Z(1);
        v(i+1) = Z(2);
    end
end