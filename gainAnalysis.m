clc; clear;
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

dataConstell = randi([0 1], 5000, M); % Случайная последовательность бит для построения сигнальных созвездий
constSNR = 60; % ОСШ для построения сигнальных созвездий

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

helper = gainAnalyserHelper();

ampl = 0.1;
ampls(:,1) = (0.05:0.05:1.5);

idealConstell = 30*ampl*qamMod(0:modOrder-1, modOrder, M);

signalIn = formSignal(constSNR, dataConstell, modOrder, h, sps, ampl, M);
rightDataIn = signalIn( (L+1)/2:end-(L+1)/2 );
samplesIn = rightDataIn(1:sps:end);

signalOut = helper.gain(signalIn);
rightDataOut = signalOut( (L+1)/2:end-(L+1)/2 );
samplesOut = rightDataOut(1:sps:end);

[evms, aclrs, paprs] = helper.calcResSaleh(constSNR, dataConstell, modOrder, h, sps, ampls, M, L, Bw, Fs);

helper.plotResSaleh(ampls, aclrs, evms, paprs, signalIn, signalOut, samplesOut, samplesIn, idealConstell, Fs);
