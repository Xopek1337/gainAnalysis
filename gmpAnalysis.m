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

% Вычисление импульсной характеристики фильтра Найквиста
h(:,1) = createH(L, Ts, T, beta);

snr = 200;
nBits = 500;

data = randi([0 1], nBits, M);

helper = gainAnalyserHelper();

ampls(:,1) = (0.2:0.2:6);
ampl = 4.8;

srrcPulse = srrcFunction(beta, sps, L);

[signalIn, modData] = formSignal(snr, data, modOrder, srrcPulse, sps, ampl, M);
signalIns = [zeros((L)/2, 1); signalIn; zeros((L)/2, 1)];

modelGMP = GMPV2();

taps = modelGMP.calcFis(signalIns, (L+1)/2);

sigOutAmplifier = taps * model_coeff;

sigOutMFilter = conv(srrcPulse, sigOutAmplifier);

sigOutMFilter = sigOutMFilter((L+1):end-(L+1));

idealConstell = ampl*qamMod(0:modOrder-1, modOrder, M);

rightDataIn = signalIn( (L+1)/2:end-(L+1)/2 );
rightDataOut = sigOutAmplifier( (L+1)/2:end-(L+1)/2 );

samplesIn = rightDataIn(1:sps:end);

samplesOut = sigOutMFilter(1:sps:end);

[evms, aclrs, paprs, phases, powers] = helper.calcResGMP(snr, data, modOrder, srrcPulse, sps, ampls, M, L, Bw, Fs, model_coeff, modelGMP);

helper.plotResGMP(aclrs, evms, paprs, rightDataIn, rightDataOut, Fs, samplesOut, idealConstell, phases, powers);
