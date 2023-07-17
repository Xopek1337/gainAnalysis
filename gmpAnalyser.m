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

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

snr = 200;
nBits = 1000;

data = randi([0 1], nBits, M);

helper = gainAnalyser();

modelGMP = GMP();

ampls(:,1) = (0.05:0.05:1.5);
ampl = 1.5;

b = [.7692, 0.1538, 0.0769, 0.0342];
a = [1,0,0,0];

signalIn = formSignal(snr, data, modOrder, h, sps, ampl, M);

signalIn2 = formSignal(snr, data, modOrder, h, sps, 0.1, M);
    
signalOut = helper.gain(signalIn);
signalOut = signalOut( (L+1)/2:end-(L+1)/2 );

[Fal, Flag, Flead] = modelGMP.calcFis(signalIn, (L+1)/2);
[Fal2, Flag2, Flead2] = modelGMP.calcFis(signalIn2, (L+1)/2);

coeffs1 = modelGMP.calcCoeffs(signalOut, Fal);
coeffs2 = modelGMP.calcCoeffs(signalOut, Flag);
coeffs3 = modelGMP.calcCoeffs(signalOut, Flead);

F = [Fal, Flag, Flead];
F2 = [Fal2, Flag2, Flead2];
coeffs = [coeffs1; coeffs2; coeffs3];

sigOut = F2 * coeffs;

[evms, aclrs, paprs, nmse] = helper.calcResGMP(snr, data, modOrder, h, sps, ampls, M, L, Bw, Fs, coeffs, modelGMP);

helper.plotResGMP(ampls, aclrs, evms, paprs, nmse, signalIn2( (L+1)/2:end-(L+1)/2 ), sigOut);