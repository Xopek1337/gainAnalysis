clc; clear;
close all

Fs = 1e7;   % частота дискретизации 
sps = 5; % число отсчетов на символ
L = 71; % длина фильтра (количество отсчетов)
T = sps/Fs;   % длительность символа
Ts = 1/Fs;  % период дискретизации
beta = 1; % степень сглаживания

M = 4; % Число бит на символ

modOrder = 2^M; % Порядок модуляции

dataConstell = randi([0 1], 5000, M); % Случайная последовательность бит для построения сигнальных созвездий
constSNR = 20; % ОСШ для построения сигнальных созвездий

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

signal = formSignal(constSNR, dataConstell, modOrder, h, sps);

b = [.7692, 0.1538, 0.0769, 0.0342];
a = [1,0,0,0];
y = filter(b, a, signal);
r = abs(y);
Ro = (30.*r)./(1+2.2.*(r.^2));
Fi = (r.^2)./(1+0.5.*(r.^2));
data_out = Ro.*exp(1j.*(Fi+angle(signal)));

rightDataOut = data_out( (L+1)/2+1:end );

samples = downsample(rightDataOut, sps);

figure(1);
plot(samples, '.');
xlabel('I'); ylabel('Q');

len=length(data_out);
n=2^nextpow2(len);
FFTY=fft(data_out,n);
FFTY=fftshift(FFTY)/len;
f1=(-n/2:n/2-1)*Fs/n;
figure(2);
plot(f1(2:length(f1)),20*log10(abs(FFTY(2:length(f1)))));
grid on;xlabel('частота, Гц');ylabel('магнитуда, дБ'); title('Спектр сигнала'); 


