result = zeros(1,N);
guess = zeros(1,N);
guess(1) = 2000;
tol = 1e-5;     %tolerância

for i=1:N-1
    val = guess(i);
    
    for j=1:N-1
        
        % resolver problema
        
    end
    
    result(i) = 1;   %resultado obtido
    
    if abs(result(i))<tol
        guess = guess(1:i);
        result = result(1:i);

        breai
    end
    
    if i>1
        m = (result(i)-result(i-1))/(guess(i)-guess(i-1));
        b = result(i)-m*guess(i);
        guess(i+1) = guess(i)-(resultado_pretendido + result(i))/m;
    end
end