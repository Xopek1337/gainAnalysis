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
nBits = 1000;

data = randi([0 1], nBits, M);

helper = gainAnalyser();

ampls(:,1) = (0.05:0.05:1.5);
ampl = 1.5;


signalIn = formSignal(snr, data, modOrder, h, sps, ampl, M);
    
modelGMP = GMPV2();

taps = modelGMP.calcFis(signalIn, (L+1)/2);

sigOut2 = taps * model_coeff;


idealConstell = ampl*qamMod(0:modOrder-1, modOrder, M);

rightDataIn = signalIn( (L+1)/2:end-(L+1)/2 );
samplesIn = rightDataIn(1:sps:end);

samplesOut = sigOut2(1:sps:end);

[evms, aclrs, paprs] = helper.calcResGMP(snr, data, modOrder, h, sps, ampls, M, L, Bw, Fs, model_coeff, v2);

helper.plotResGMP(ampls, aclrs, evms, paprs, circshift(rightDataIn,1), sigOut2, Fs, samplesOut, samplesIn, idealConstell);
