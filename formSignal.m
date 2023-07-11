function [signal] = formSignal(snr, data, modOrder, h, sps, ampl, M)
    dataSym = bi2de(data);
    modData = qamMod(dataSym, modOrder, M);
    
    output = zeros(length(data)*sps, 1);
    output(1) = modData(1);
    for i = 1:length(modData)-1
       output(sps*i+1) = modData(i+1); 
    end
    
    signal = ampl * conv(h, output);
    varSignal = var(signal);

    varNoise = varSignal*10^(-snr/10);
    signal = signal+sqrt(varNoise/2)*(randn(size(signal))+1i*randn(size(signal)));
end
