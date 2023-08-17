function [p] = srrcFunction(beta, sps, L)
    Tsym = sps; t=-(L/2):1:(L/2);

    num = sin(pi*t*(1-beta)/Tsym)+((4*beta*t/Tsym).*cos(pi*t*(1+beta)/Tsym));
    den = pi*t.*(1-(4*beta*t/Tsym).^2)/Tsym;
    p = 1/sqrt(Tsym)*num./den;

    p(ceil(length(p)/2))=1/sqrt(Tsym)*((1-beta)+4*beta/pi);

    temp=(beta/sqrt(2*Tsym))*( (1+2/pi)*sin(pi/(4*beta)) + (1-2/pi)*cos(pi/(4*beta)));
    
    p(t==Tsym/(4*beta))=temp; p(t==-Tsym/(4*beta))=temp;
end
