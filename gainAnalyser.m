classdef gainAnalyser
    properties
        a {mustBeNumeric}
        b {mustBeNumeric}
    end
    methods
        function obj = gainAnalyser(a, b)
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
            aclr = 10*log10(sum(spm((Fs/2 - Bw/2)/Fs*n : (Fs/2 + Bw/2)/Fs*n)) / sum(spm((Fs/2 - 3*Bw/2)/Fs*n : (Fs/2 - Bw/2)/Fs*n)));
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
            xlabel('Амплитуда');ylabel('ACLR, дБ');
            
            figure(2);
            plot(ampls, evms);
            xlabel('Амплитуда');ylabel('EVM, дБ');
            
            figure(3);
            plot(ampls, paprs);
            xlabel('Амплитуда на входе');ylabel('PAPR, дБ');
            
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

        function [evms, aclrs, paprs] = calcResGMP(obj, constSNR, dataConstell, modOrder, h, sps, ampls, M, L, Bw, Fs, coeffs, modelGMP)
            for i=1:length(ampls)
                sig = formSignal(constSNR, dataConstell, modOrder, h, sps, ampls(i), M);
                
                idealCons = ampls(i)*qamMod(0:15, modOrder, M);
                
                F = modelGMP.calcFis(sig, (L+1)/2);
                sigOutGMP = F * coeffs;

                sigOutSaleh = obj.gain(sig);
                sigOutSaleh = sigOutSaleh( (L+1)/2:end-(L+1)/2 );

                len=length(sigOutGMP);
                n=2^nextpow2(len);
                
                samplesOut = sigOutGMP(1:sps:end);
                
                FFTYOut = obj.calcSpectrum(sigOutGMP);
                
                evms(i) = obj.calcEVM(samplesOut, idealCons);
                aclrs(i) = obj.calcACLR(abs(FFTYOut).^2, n, Bw, Fs);
                paprs(i) = obj.calcPAPR(sigOutGMP);
            end
        end

        function [] = plotResGMP(obj, ampls, aclrs, evms, paprs, sigIn, sigOut, Fs, samplesOut, samplesIn, idealConstell)
            figure(1);
            plot(ampls, aclrs);
            xlabel('Амплитуда');ylabel('ACLR, дБ');
            
            figure(2);
            plot(ampls, evms);
            xlabel('Амплитуда');ylabel('EVM, дБ');
            
            figure(3);
            plot(ampls, paprs);
            xlabel('Амплитуда на входе');ylabel('PAPR, дБ');
            
            figure(4);
            semilogy((0:0.05:15), obj.calcCCDF(sigOut));
            xlabel('PAPR, дБ');ylabel('CCDF');

            figure(5);
            plot(abs(sigIn), abs(sigOut), '.');
            xlabel('Амплитуда на входе');ylabel('Амплитуда на выходе GMP усилителя');
            
            phase = angle(sigOut) - angle(sigIn);
            phase = sin(phase);
            phase = asin(phase);
            %invPhase = find(phase<-6.1);
            %phase(invPhase) = phase(invPhase) + 2*pi;
            
            figure(6);
            plot(abs(sigIn), phase, '.');
            xlabel('Амплитуда на входе');ylabel('Фаза на выходе GMP усилителя');

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
            hold on;
            plot(idealConstell,'*');
            plot(samplesIn, '.');
            hold off;
            xlabel('I'); ylabel('Q'); legend('Символы на выходе усилителя', 'Символы на идеальном созвездии', 'Символы при линейном усилении');
        end
    end
end
