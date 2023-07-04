function [signal] = formSignal(snr, data, modOrder, h, sps)
    dataSym = bi2de(data);
    modData = qammod(dataSym, modOrder, 'UnitAveragePower' , true);

    output = upsample(modData, sps);
    
    signal = conv(h, output);
    varSignal = var(signal);

    varNoise = varSignal*10^(-snr/10);
    signal = signal+sqrt(varNoise/2)*(randn(size(signal))+1i*randn(size(signal)));
end