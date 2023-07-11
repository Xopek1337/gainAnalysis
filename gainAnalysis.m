clc; clear;
close all

Fs = 1e7;   % частота дискретизации 
sps = 10; % число отсчетов на символ
L = 71; % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 1; % степень сглаживания
Bw = (1+beta) / T; % ширина главного лепестка

M = 4; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

dataConstell = randi([0 1], 5000, M); % Случайная последовательность бит для построения сигнальных созвездий
constSNR = 60; % ОСШ для построения сигнальных созвездий

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

ampl = 0.1;
ampls(:,1) = (0.05:0.05:1.5);

ros = [];
fis = [];
evms = [];
aclrs = [];
paprs = [];
ccdf = [];

idealConstell = 30*ampl*qamMod(0:modOrder-1, modOrder, M);

signal = formSignal(constSNR, dataConstell, modOrder, h, sps, ampl, M);

b = [.7692, 0.1538, 0.0769, 0.0342];
a = [1,0,0,0];

[data_out, r] = gain(a, b, signal);

rightDataOut = data_out( (L+1)/2:end-(L+1)/2 );

samples = rightDataOut(1:sps:end);

rightDataOut2 = signal( (L+1)/2:end-(L+1)/2 );

samples2 = rightDataOut2(1:sps:end);

len=length(data_out);
n=2^nextpow2(len);

FFTY = calcSpectrum(data_out);
FFTY2 = calcSpectrum(signal);

f1=(-n/2:n/2-1)*Fs/n;
figure(1);
FFTY = medfilt1(abs(FFTY),150,'truncate');
FFTY2 = medfilt1(abs(FFTY2),150,'truncate');

plot(f1(2:length(f1)), 20*log10(abs(FFTY(2:length(f1)))));
hold on;
plot(f1(2:length(f1)), 20*log10(abs(FFTY2(2:length(f1))))+29.5);
hold off;
grid on; xlabel('частота, Гц'); ylabel('магнитуда, дБ'); title('Спектр сигнала'); legend('Спектр на выходе усилителя','Спектр на входе усилителя');

figure(2);
plot(samples, '.');

hold on;
plot(idealConstell,'*');
plot(samples2*30, '.');
hold off;
xlabel('I'); ylabel('Q'); legend('Символы на выходе усилителя', 'Символы на идеальном созвездии', 'Символы при линейном усилении');

for i=1:length(ampls)
    sig = formSignal(constSNR, dataConstell, modOrder, h, sps, ampls(i), M);
    idealCons = 30*ampls(i)*qamMod(0:15, modOrder, M);
    
    [sigOut, ro, fi] = gain(a, b, sig);
    
    rightSigOut = sigOut( (L+1)/2:end-(L+1)/2 );

    samplesOut = rightSigOut(1:sps:end);
    
    FFTYOut = calcSpectrum(sigOut);
    
    evms(i) = calcEVM(samplesOut, idealCons);
    aclrs(i) = calcACLR(abs(FFTYOut).^2, n, Bw, Fs);
    ros(i) = sum(ro)/length(ro);
    fis(i) = sum(fi)/length(fi);
    paprs(i) = calcPAPR(sigOut);    
end

figure(3);
plot(ampls, aclrs);
xlabel('Амплитуда');ylabel('ACLR, дБ');

figure(4);
plot(ampls, evms);
xlabel('Амплитуда');ylabel('EVM, дБ');

figure(5);
plot(ampls, fis);
xlabel('Амплитуда на входе');ylabel('Фаза на выходе');

figure(6);
plot(ampls, ros);
xlabel('Амплитуда на входе');ylabel('Амплитуда на выходе');

figure(7);
plot(ampls, paprs);
xlabel('Амплитуда на входе');ylabel('PAPR, дБ');

figure(8);
semilogy((0:0.05:15), calcCCDF(data_out));
xlabel('PAPR, дБ');ylabel('CCDF');

function [papr, ccdf] = calcPAPR(sig)
    P = abs(sig).^2;
    papr = 10*log10((max(P)/mean(P)));
end

function [ccdf] = calcCCDF(sig)
    ccdf = [];
    papr = (0:0.05:15);
    P = abs(sig).^2;
    Pr = 10*log10(P/mean(P));
    for i = 1:length(papr)
        ccdf(i) = length(find(Pr >= papr(i)))/length(sig);
    end
end

function [aclr] = calcACLR(spm, n, Bw, Fs)
    aclr = 10*log10(sum(spm((Fs/2 - Bw/2)/Fs*n : (Fs/2 + Bw/2)/Fs*n)) / sum(spm((Fs/2 - 3*Bw/2)/Fs*n : (Fs/2 - Bw/2)/Fs*n)));
end

function [evm] = calcEVM(samples, constell)
    sumCnst = 0;
    evm = [];
    for i=1:length(samples)
       [val, idx]=min(samples(i) - constell);
       evm(i) = abs(val^2);
    end
    
    for i=1:length(constell)
       sumCnst = sumCnst + abs(constell(i))^2;
    end
    
    sumCnst = sumCnst/length(constell);
    evm = 20*log10(sqrt((sum(evm) / length(samples))/sumCnst));
end

function [spm] = calcSpectrum(data_out)
    len=length(data_out);
    n=2^nextpow2(len);
    FFTY=fft(data_out,n);
    spm=fftshift(FFTY)/len;
end

function [sig, Ro, Fi] = gain(a, b, signal)
    y = filter(b, a, signal);
    r = abs(y);
    Ro = (30.*r)./(1+2.2.*(r.^2));
    Fi = (r.^2)./(1+0.5.*(r.^2));
    
    sig = Ro.*exp(1j.*(Fi+angle(signal)));
end
