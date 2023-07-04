function [h] = createH(L, Ts, T, beta) % фильтр Найквиста типа 'Приподнятый косинус'
    h = zeros(1, L); 
    if mod(L, 2)==0     
        M = L/2;
    else
        M = (L-1)/2;
    end
    
    for n = -M:M     
        num = sin(pi*n*Ts/T)*cos(beta*pi*n*Ts/T);     
        den = (pi*n*Ts/T)*(1-(2*beta*n*Ts/T)^2);     
        h(n+M+1) = num/den;     
        if (1-(2*beta*n*Ts/T)^2) == 0         
            h(n+M+1) = pi/4*sin(pi*n*Ts/T)/(pi*n*Ts/T);     
        end
        if n == 0         
            h(n+M+1) = cos(beta*pi*n*Ts/T)/(1-(2*beta*n*Ts/T)^2); % sinc == 1    
        end
    end
end