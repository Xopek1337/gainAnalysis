classdef GMPV2
    properties
        order {mustBeNumeric}
        i0_max {mustBeNumeric}
        i1_max {mustBeNumeric}
        depth {mustBeNumeric}
    end
    methods
        function obj = GMPV2(order, i0_max, i1_max)
            if nargin == 3
                obj.order = order;
                obj.i0_max = i0_max;
                obj.i1_max = i1_max;
            else
                obj.order = 5;
                obj.i0_max = 11;
                obj.i1_max = 11;
            end
            obj.depth = max(obj.i0_max, obj.i1_max);
        end

        function poly_coeffs = calcCoeffs(obj, y, Fis)
            poly_coeffs = inv(Fis'*Fis + 0.1 * eye(size(Fis'*Fis)))*Fis'*y;
        end

        function Fal = calcFis(obj, x, start)
            Fal = zeros(length(x)-2*start, obj.depth + (obj.order-1)*obj.i0_max*obj.i1_max);

            for i = start:length(x)-start
                Fal(i-start+1, :)= featureGenerationGMPmodel(obj, x, i, 1);
            end
        end

        function taps = featureGenerationGMPmodel(obj, signal, n, allPowers)
            if n - obj.depth < 0
                tapsDelay = [zeros(1, abs(n-obj.depth)), signal(1:n)];
            elseif n > length(signal)
                tapsDelay = [signal(n-obj.depth+1:length(signal)), zeros(1, abs(length(signal)-n))];
            else
                tapsDelay = signal(n-obj.depth+1:n);
            end
            
            if allPowers
                taps = zeros(obj.depth + (obj.order-1)*obj.i0_max*obj.i1_max, 1);
                orderStep = 1;
            else
                taps = zeros(obj.depth + floor((obj.order-1)/2)*obj.i0_max*obj.i1_max, 1);
                orderStep = 2;
            end
            
            j = 1;
            for order = 1 : orderStep : obj.order
                if order == 1
                    for i0 = 1 : obj.i0_max
                        taps(j) = tapsDelay(length(tapsDelay)-i0+1);
                        j = j + 1;
                    end
                else
                    for i0 = 1 : obj.i0_max
                        for i1 = 1 : obj.i1_max
                            taps(j) = tapsDelay(length(tapsDelay)-i0+1)*(abs(tapsDelay(length(tapsDelay)-i1+1)) ^ (order-1));
                            j = j + 1;
                        end
                    end
                end
            end
        end
    end
end
