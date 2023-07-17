classdef GMP
    properties
        order {mustBeNumeric}
        memory_depth {mustBeNumeric}
        P {mustBeNumeric}
        Q {mustBeNumeric}
    end
    methods
        function obj = GMP(memory_depth, order, Q, P)
            if nargin == 4
                obj.memory_depth = memory_depth;
                obj.order = order;
                obj.P = Q;
                obj.Q = P;
            else
                obj.memory_depth = 7;
                obj.order = 4;
                obj.P = 0;
                obj.Q = 0;
            end
        end

        function poly_coeffs = calcCoeffs(obj, y, Fis)
            poly_coeffs = inv(Fis'*Fis)*Fis'*y;
        end

        function [Fal, Flag, Flead] = calcFis(obj, x, start)
            Fal = zeros(length(x)-2*start, obj.memory_depth*obj.order);
            Flag = zeros(length(x)-2*start, obj.memory_depth*obj.order*(obj.P));
            Flead = zeros(length(x)-2*start, obj.memory_depth*obj.order*(obj.Q));

            for i = start:length(x)-start
                [Fal(i-start+1,:), Flag(i-start+1,:), Flead(i-start+1,:)] = setFi(obj, x, i);
            end
        end

        function [Falign, Flag, Flead] = setFi(obj, x, ind)
            Falign(:,1) = zeros(obj.memory_depth*obj.order, 1);
            Flag(:,1) = zeros(obj.memory_depth*obj.order*obj.P, 1);
            Flead(:,1) = zeros(obj.memory_depth*obj.order*obj.Q, 1);

            for i = 1:obj.memory_depth
                for j = 1:obj.order
                    Falign(j + (i-1)*obj.order) = x(ind + 1 - i) * abs(x(ind + 1 - i))^(j-1);
                    for k = 1:obj.P
                        Flag(k + (j-1)*obj.P + (i-1)*obj.order*obj.P) = x(ind + 1 - i) * abs(x(ind - i - k - 1))^(j-1);
                    end
                    for k = 1:obj.Q
                        Flead(k + (j-1)*obj.Q + (i-1)*obj.order*obj.P) = x(ind + 1 - i) * abs(x(ind - i + k - 1))^(j-1);
                    end
                end
            end

            Falign = Falign.';
            Flag = Flag.';
            Flead = Flead.';
        end
    end
end