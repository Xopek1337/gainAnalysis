%clc; clear;
close all

Fs = 1e7;   % частота дискретизации 
sps = 20; % число отсчетов на символ
L = 72; % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 1; % степень сглаживания
Bw = (1+beta) / T; % ширина главного лепестка

M = 4; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

snr = 200;
nBits = 500;

data = randi([0 1], nBits, M);

helper = gainAnalyserHelper();

ampl = 5.1;

srrcPulse = srrcFunction(beta, sps, L);

signal = formSignal(snr, data, modOrder, srrcPulse, sps, ampl, M);

power = (sum(abs(signal .^2)) / length(signal));

modelGMP = GMPV2();

signals = [zeros((L)/2, 1); signal; zeros((L)/2, 1)];

taps = modelGMP.calcFis(signals, (L+1)/2);

sigOutAmplifier = taps * model_coeff;

sigOutFiltered = conv(srrcPulse, sigOutAmplifier);

sigOutFiltered = sigOutFiltered((L+1):end-(L+1));

idealConstell = qamMod(0:15, 16, 4);

start = pi/4;
phShift = start;
evm = 0;

samplesOut = sigOutFiltered(1:sps:end);

for i = 1:8000
    phases(i) = phShift;

    sigOutShifted = abs(sigOutFiltered).*exp(1j.*(angle(sigOutFiltered) - phShift));

    phShift = start-i/5000;
    
    samplesOutShifted = sigOutShifted(1:sps:end);
    
    samplesOutShifted = samplesOutShifted/sqrt(sum(abs(samplesOutShifted .^2)) / length(samplesOutShifted));

    evms(i) = calcEVM(samplesOutShifted, idealConstell);
end

[valMin, ind] = min(evms);

phShift = start-ind/5000;

sigOutShifted = abs(sigOutFiltered).*exp(1j.*(angle(sigOutFiltered) - phShift));
    
samplesOutShifted = sigOutShifted(1:sps:end);

samplesOutShifted = samplesOutShifted/sqrt(sum(abs(samplesOutShifted .^2)) / length(samplesOutShifted));

samplesOut = samplesOut/sqrt(sum(abs(samplesOut .^2)) / length(samplesOut));

figure(1);
plot(samplesOut, '.');
grid on;
hold on;
plot(samplesOutShifted, '.');
plot(idealConstell, '*'); legend('Символы на выходе усилителя без сдвига на среднюю фазу', 'Символы на выходе усилителя со сдвигом на среднюю фазу', 'Символы на идеальном созвездии');
hold off;

figure(2);
plot(phases, evms);
grid on;
xlabel('Фазовый сдвиг, рад'); ylabel('EVM, дБ');


rightDataIn = signal( (L+1)/2:end-(L+1)/2 );
rightDataOut = sigOutAmplifier( (L+1)/2:end-(L+1)/2 );

figure(3);
plot(abs(rightDataIn), abs(rightDataOut), '.');
grid on;
xlabel('Амплитуда на входе');ylabel('Амплитуда на выходе GMP усилителя');

phase = angle(rightDataOut) - angle(rightDataIn);
phase = sin(phase);
phase = asin(phase);

figure(4);
plot(abs(rightDataIn), phase, '.');
grid on;
xlabel('Амплитуда на входе');ylabel('Фаза на выходе GMP усилителя');

evm2 = calcEVM(samplesOutShifted, idealConstell);

 len=length(rightDataOut);
            n=2^nextpow2(len);
            
            FFTYOut = helper.calcSpectrum(rightDataOut);
            FFTYIn = helper.calcSpectrum(rightDataIn);
            
            f1=(-n/2:n/2-1)*Fs/n;

            figure(7);
            FFTYOut = medfilt1(abs(FFTYOut),150,'truncate');
            FFTYIn = medfilt1(abs(FFTYIn),150,'truncate');
            
            plot(f1(2:length(f1)), 20*log10(abs(FFTYOut(2:length(f1)))));
            hold on;
            plot(f1(2:length(f1)), 20*log10(abs(FFTYIn(2:length(f1)))));
            hold off;
            grid on; xlabel('частота, Гц'); ylabel('магнитуда, дБ'); title('Спектр сигнала'); legend('Спектр на выходе усилителя','Спектр на входе усилителя');

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
