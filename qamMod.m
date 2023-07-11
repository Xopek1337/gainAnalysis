function symbols = qamMod(x, M, b)
    b_I = floor(b/2);
    M_I = 2^b_I;
    b_Q = ceil(b/2);
    M_Q = 2^b_Q;
    x_I = -(M_I-1)/2:(M_I-1)/2;
    x_Q = -(M_Q-1)/2:(M_Q-1)/2;
    x_IQ = x_I + 1i*x_Q.';
    x_IQ = x_IQ / sqrt((M_I^2+M_Q^2-2)/12);

    symbols = x_IQ(x+1);   
end
