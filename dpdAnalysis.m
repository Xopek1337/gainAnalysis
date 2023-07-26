%clc; clear;
close all

Fs = 1e7;   % частота дискретизации 
sps = 20; % число отсчетов на символ
L = 71; % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 1; % степень сглаживания
Bw = (1+beta) / T; % ширина главного лепестка

M = 4; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

snr = 200;
nBits = 500;

data = randi([0 1], nBits, M);

helper = gainAnalyser();

ampl = 1.2;

b = [.7692, 0.1538, 0.0769, 0.0342];
a = [1,0,0,0];

signalIn = formSignal(snr, data, modOrder, h, sps, ampl, M);

v2 = GMPV2();

taps2 = v2.calcFis(signalIn, (L+1)/2);

sigOut2 = (taps2 * model_coeff);

idealConstell = ampl*qamMod(0:modOrder-1, modOrder, M);

rightDataIn = signalIn( (L+1)/2:end-(L+1)/2 );

samplesOut = sigOut2(1:sps:end);

dpd = GMPV2(5, 5, 5);

taps3 = dpd.calcFis(signalIn, (L+1)/2);

sigOutPAs = [zeros((L+1)/2 - 1, 1); sigOut2; zeros((L+1)/2, 1)];
tapsPA = dpd.calcFis(sigOutPAs, (L+1)/2);
coeffss = dpd.calcCoeffs(rightDataIn, tapsPA);
sigOutDPD = taps3 * coeffss;
sigOutDPDs = [zeros((L+1)/2 - 1, 1); sigOutDPD; zeros((L+1)/2, 1)];
tapsDPD = v2.calcFis(sigOutDPDs, (L+1)/2);
sigOutPA = tapsDPD * model_coeff;

samplesIn = sigOutPA(1:sps:end);
samplesIn2 = rightDataIn(1:sps:end);
evm = helper.calcEVM(samplesIn, idealConstell);
disp(evm);

plot(samplesIn, '.');
hold on;
plot(samplesOut, '.');
samplesIn2 = rightDataIn(1:sps:end);
plot(samplesIn2, '*'); legend('Символы на выходе усилителя c dpd', 'Символы на выходе усилителя без dpd', 'Символы на входе');
hold off;

len=length(sigOutPA);
n=2^nextpow2(len);
            
FFTYOut = helper.calcSpectrum(sigOutPA);
FFTYIn = helper.calcSpectrum(sigOut2);
FFTYIn2 = helper.calcSpectrum(rightDataIn);
            
f1=(-n/2:n/2-1)*Fs/n;

figure(2);
FFTYOut = medfilt1(abs(FFTYOut),150,'truncate');
FFTYIn = medfilt1(abs(FFTYIn),150,'truncate');
FFTYIn2 = medfilt1(abs(FFTYIn2),150,'truncate');
            
plot(f1(2:length(f1)), 20*log10(abs(FFTYOut(2:length(f1)))));
hold on;
plot(f1(2:length(f1)), 20*log10(abs(FFTYIn(2:length(f1)))));
plot(f1(2:length(f1)), 20*log10(abs(FFTYIn2(2:length(f1)))));
hold off;
grid on; xlabel('частота, Гц'); ylabel('магнитуда, дБ'); title('Спектр сигнала'); legend('Спектр на выходе усилителя c DPD','Спектр на выходе усилителя без DPD', 'Спектр на входе усилителя');

figure(3);
plot(abs(rightDataIn), abs(sigOutPA), '.');
xlabel('Амплитуда на входе с DPD'); ylabel('Амплитуда на выходе с DPD');

figure(4);
plot(abs(rightDataIn), abs(sigOut2), '.');
xlabel('Амплитуда на входе без DPD'); ylabel('Амплитуда на выходе без DPD');