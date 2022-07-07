classdef multipliers < handle
    properties % (Access=private)
        cache
        shape
        res
    end
    methods(Static) % (Access=private)
        function out=linspace_near_0(len, pos_0, edge)
            % Give 1d arr of linespace [-0.5, 0.5)
            % Example of edge: 
            % true->0, 0.25, -0.5, 0
            % false->0.125, 0.375, -0.375, -0.125
            arguments
                len, pos_0=0, edge=true
            end
            out = 0:len-1;
            out = mod((out + 0.5*~edge)/len - pos_0 + 0.5, 1) - 0.5;
        end
    end

    methods
        function obj=multipliers(shape, res)
            obj.shape = shape;
            obj.res = res;
        end
        function arr=gamma(obj)
            if ~isfield(obj.cache, 'gamma')
                EPS = 1e-10;
                alpha = obj.linspace_near_0(obj.shape(2)) / obj.res(1);
                beta = obj.linspace_near_0(obj.shape(1)) / obj.res(2);
                obj.cache.gamma = sqrt(max(1-alpha.^2-beta'.^2, EPS));
            end
            arr = obj.cache.gamma;
        end
        function [arr, alpha_beta]=tilt(obj, alpha_beta, trunc)
            arguments 
                obj, alpha_beta, trunc=true
            end
            nm = obj.shape([2 1]) .* obj.res([1 2]);
            if trunc
                alpha_beta = fix(alpha_beta .* nm) ./ nm;
            end
%             if ~isfield(obj.cache, 'tilt') && obj.cache.tilt(alpha_beta)
%             end
            xr = linspace(0,1,obj.shape(2)+1) * alpha_beta(1) * nm(1);
            yr = linspace(0,1,obj.shape(1)+1) * alpha_beta(2) * nm(2);
            xr = xr(1:end-1);
            yr = yr(1:end-1);
            phase = mod(xr - yr', 1);
            arr = exp(2j * pi * phase);
            % normalize the phase using center point
            arr = arr / arr(floor(obj.shape(1)/2), floor(obj.shape(2)/2));
        end
        function arr=soft_crop(obj, width, strength, total_slices)
            fac1 = log(100 * strength) + 0.8;
            fac2 = 100 * strength / total_slices;
            x_fac = exp(-exp(-(obj.linspace_near_0(obj.shape(2), 0) * 2 / width) .^ 2 * fac1)*fac2);
            y_fac = exp(-exp(-(obj.linspace_near_0(obj.shape(1), 0) * 2 / width) .^ 2 * fac1)*fac2);
            arr = x_fac .* y_fac';
        end
        function arr=gaussian(obj, mu, sigma)
            % mu in 0~1, sigma in px
            x_fac = exp(-(obj.linspace_near_0(obj.shape(2), mu(1))) .^ 2 / 2 / (sigma(1) / obj.shape(2)) ^ 2);
            y_fac = exp(-(obj.linspace_near_0(obj.shape(1), mu(2))) .^ 2 / 2 / (sigma(2) / obj.shape(1)) ^ 2);
            arr = x_fac .* y_fac';
        end
    end
end