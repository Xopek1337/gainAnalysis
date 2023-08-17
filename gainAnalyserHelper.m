classdef gainAnalyserHelper
    properties
        a {mustBeNumeric}
        b {mustBeNumeric}
    end
    methods
        function obj = gainAnalyserHelper(a, b)
            if nargin == 2
                obj.a = a;
                obj.b = b;
            else
                obj.a = [1,0,0,0];
                obj.b = [.7692, 0.1538, 0.0769, 0.0342];
            end
        end

        function [papr, ccdf] = calcPAPR(obj, sig)
            P = abs(sig).^2;
            papr = 10*log10((max(P)/mean(P)));
        end
        
        function [ccdf] = calcCCDF(obj, sig)
            ccdf = [];
            papr = (0:0.05:15);
            P = abs(sig).^2;
            Pr = 10*log10(P/mean(P));
            for i = 1:length(papr)
                ccdf(i) = length(find(Pr >= papr(i)))/length(sig);
            end
        end
        
        function [aclr] = calcACLR(obj, spm, n, Bw, Fs)
            aclr = 10*log10(sum(spm((Fs/2 - Bw/2)/Fs*n : (Fs/2 + Bw/2)/Fs*n)) / (sum(spm((Fs/2 - 3*Bw/2)/Fs*n : (Fs/2 - Bw/2)/Fs*n)) + sum(spm((Fs/2 + Bw/2)/Fs*n : (Fs/2 + 3*Bw/2)/Fs*n))));
        end
        
        function [evm] = calcEVM(obj, samples, constell)
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
        
        function [spm] = calcSpectrum(obj, data_out)
            len=length(data_out);
            n=2^nextpow2(len);
            FFTY=fft(data_out,n);
            spm=fftshift(FFTY)/len;
        end
        
        function [sig, Ro, Fi] = gain(obj, signal)
            y = filter(obj.b, obj.a, signal);
            r = abs(y);
            Ro = (30.*r)./(1+2.2.*(r.^2));
            Fi = (r.^2)./(1+0.5.*(r.^2));
            
            sig = Ro.*exp(1j.*(Fi+angle(signal)));
        end

        function [evms, aclrs, paprs] = calcResSaleh(obj, constSNR, dataConstell, modOrder, h, sps, ampls, M, L, Bw, Fs)
            for i=1:length(ampls)
                sig = formSignal(constSNR, dataConstell, modOrder, h, sps, ampls(i), M);
                
                idealCons = 30*ampls(i)*qamMod(0:15, modOrder, M);
                
                sigOut = obj.gain(sig);

                len=length(sigOut);
                n=2^nextpow2(len);
                
                rightSigOut = sigOut( (L+1)/2:end-(L+1)/2 );
            
                samplesOut = rightSigOut(1:sps:end);
                
                FFTYOut = obj.calcSpectrum(sigOut);
                
                evms(i) = obj.calcEVM(samplesOut, idealCons);
                aclrs(i) = obj.calcACLR(abs(FFTYOut).^2, n, Bw, Fs);
                paprs(i) = obj.calcPAPR(sigOut);    
            end
        end

        function [] = plotResSaleh(obj, ampls, aclrs, evms, paprs, sigIn, sigOut, samplesOut, samplesIn, idealConstell, Fs)
            figure(1);
            plot(ampls, aclrs);
            xlabel('Средняя мощность входного сигнала');ylabel('ACLR, дБ');
            
            figure(2);
            plot(ampls, evms);
            xlabel('Средняя мощность входного сигнала');ylabel('EVM, дБ');
            
            figure(3);
            plot(ampls, paprs);
            xlabel('Средняя мощность входного сигнала');ylabel('PAPR, дБ');
            
            figure(4);
            semilogy((0:0.05:15), obj.calcCCDF(sigOut));
            xlabel('PAPR, дБ');ylabel('CCDF');

            len=length(sigOut);
            n=2^nextpow2(len);
            
            FFTYOut = obj.calcSpectrum(sigOut);
            FFTYIn = obj.calcSpectrum(sigIn);
            
            f1=(-n/2:n/2-1)*Fs/n;

            figure(5);
            FFTYOut = medfilt1(abs(FFTYOut),150,'truncate');
            FFTYIn = medfilt1(abs(FFTYIn),150,'truncate');
            
            plot(f1(2:length(f1)), 20*log10(abs(FFTYOut(2:length(f1)))));
            hold on;
            plot(f1(2:length(f1)), 20*log10(abs(FFTYIn(2:length(f1))))+29.5);
            hold off;
            grid on; xlabel('частота, Гц'); ylabel('магнитуда, дБ'); title('Спектр сигнала'); legend('Спектр на выходе усилителя','Спектр на входе усилителя');

            figure(6);
            plot(samplesOut, '.');
            hold on;
            plot(idealConstell,'*');
            plot(samplesIn*30, '.');
            hold off;
            xlabel('I'); ylabel('Q'); legend('Символы на выходе усилителя', 'Символы на идеальном созвездии', 'Символы при линейном усилении');
            
            figure(7);
            plot(abs(sigIn), abs(sigOut), '.');
            xlabel('Амплитуда на входе');ylabel('Амплитуда на выходе GMP усилителя');
            
            phase = angle(sigOut) - angle(sigIn);
            invPhase = find(phase<-2*pi+max(phase));
            phase(invPhase) = phase(invPhase) + 2*pi;
            
            figure(8);
            plot(abs(sigIn), phase, '.');
            xlabel('Амплитуда на входе');ylabel('Фаза на выходе GMP усилителя');
        end

        function [evms, aclrs, paprs, phases, powers] = calcResGMP(obj, constSNR, dataConstell, modOrder, h, sps, ampls, M, L, Bw, Fs, coeffs, modelGMP)
            for i=1:length(ampls)
                [sig, modData] = formSignal(constSNR, dataConstell, modOrder, h, sps, ampls(i), M);
                power =  (sum(abs(sig .^2)) / length(sig));

                sig = [zeros((L)/2, 1); sig; zeros((L)/2, 1)];
                
                idealCons = qamMod(0:15, modOrder, M);
                
                F = modelGMP.calcFis(sig, (L+1)/2);
                sigOutGMP = F * coeffs;

                sigOutGMPFiltered = conv(h, sigOutGMP);

                sigOutGMPFiltered = sigOutGMPFiltered((L+1):end-(L+1));

                len=length(sigOutGMP);
                n=2^nextpow2(len);
                
                samplesOut = sigOutGMPFiltered(1:sps:end);

                start = pi/4;

                phDelay = start;

                for j = 1:8000
                    sigOutShifted = abs(sigOutGMPFiltered).*exp(1j.*(angle(sigOutGMPFiltered) - phDelay));

                    phDelay = start-j/5000;
                    
                    samplesOutShifted = sigOutShifted(1:sps:end);
              
                    samplesOutShifted = samplesOutShifted/sqrt(sum(abs(samplesOutShifted .^2)) / length(samplesOutShifted));
                
                    evmsShifted(j) = obj.calcEVM(samplesOutShifted, idealCons);
                end
                [valMin, ind] = min(evmsShifted);

                phDelay = start-ind/5000;
                
                sigOutShifted = abs(sigOutGMPFiltered).*exp(1j.*(angle(sigOutGMPFiltered) - phDelay));

                samplesOutShifted = sigOutShifted(1:sps:end);

                samplesOutShifted = samplesOutShifted/sqrt(sum(abs(samplesOutShifted .^2)) / length(samplesOutShifted));
                
                FFTYOut = obj.calcSpectrum(sigOutGMP);

                powers(i) = power;
                phases(i) = phDelay;
                evms(i) = obj.calcEVM(samplesOutShifted, idealCons);
                aclrs(i) = obj.calcACLR(abs(FFTYOut).^2, n, Bw, Fs);
                paprs(i) = obj.calcPAPR(sigOutGMP);
            end
        end

        function [] = plotResGMP(obj, aclrs, evms, paprs, sigIn, sigOut, Fs, samplesOut, idealConstell, phases, powers)
            figure(1);
            plot(powers, aclrs);
            grid on;
            xlabel('Средняя мощность входного сигнала');ylabel('ACLR, дБ');
            
            figure(2);
            plot(powers, evms);
            grid on;
            xlabel('Средняя мощность входного сигнала');ylabel('EVM, дБ');
            
            figure(3);
            plot(powers, paprs);
            grid on;
            xlabel('Средняя мощность входного сигнала');ylabel('PAPR, дБ');
            
            figure(4);
            semilogy((0:0.05:15), obj.calcCCDF(sigOut));
            grid on;
            xlabel('дБ выше средней мощности');ylabel('CCDF');

            figure(5);
            plot(abs(sigIn), abs(sigOut), '.');
            grid on;
            xlabel('Мгновенная амплитуда на входе');ylabel('Мгновенная амплитуда на выходе GMP усилителя');
            
            phase = angle(sigOut) - angle(sigIn);
            phase = sin(phase);
            phase = asin(phase);
            %invPhase = find(phase<-6.1);
            %phase(invPhase) = phase(invPhase) + 2*pi;
            
            figure(6);
            plot(abs(sigIn), phase, '.');
            grid on;
            xlabel('Мгновенная амплитуда на входе');ylabel('Мгновенная фаза на выходе GMP усилителя');

            len=length(sigOut);
            n=2^nextpow2(len);
            
            FFTYOut = obj.calcSpectrum(sigOut);
            FFTYIn = obj.calcSpectrum(sigIn);
            
            f1=(-n/2:n/2-1)*Fs/n;

            figure(7);
            FFTYOut = medfilt1(abs(FFTYOut),150,'truncate');
            FFTYIn = medfilt1(abs(FFTYIn),150,'truncate');
            
            plot(f1(2:length(f1)), 20*log10(abs(FFTYOut(2:length(f1)))));
            hold on;
            plot(f1(2:length(f1)), 20*log10(abs(FFTYIn(2:length(f1)))));
            hold off;
            grid on; xlabel('частота, Гц'); ylabel('магнитуда, дБ'); title('Спектр сигнала'); legend('Спектр на выходе усилителя','Спектр на входе усилителя');

            figure(8);
            plot(samplesOut, '.');
            grid on;
            hold on;
            plot(idealConstell,'*');
            hold off;
            xlabel('I'); ylabel('Q'); legend('Символы на выходе усилителя', 'Символы на идеальном созвездии');
        end
    end
end
