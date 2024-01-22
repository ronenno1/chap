classdef ffts   
    methods(Static)
        function data = keap_only(data, order)
            center = floor(length(data)/2+1);
            N_clear = center-order-1;
            S_cleared = fftshift(fft(data));   
            S_cleared(1:N_clear) = 0;
            S_cleared(end-N_clear+2:end) = 0;
            S_cleared = fftshift(S_cleared);
            data = ifft(S_cleared);  
        end
        
        function [P1, f] = do_fft(data, Fs)
            norm_data = data - mean(data, 'omitnan');
            norm_data(isnan(norm_data)) = 0;
            for signal=1:size(norm_data, 1)
                Y(signal, :) = fft(norm_data(signal, :))-1;
            end
            L = length(data);
            P2 = abs(Y/L);
            P1 = P2(:, 1:floor(L/2)+1);
            P1(:, 2:end-1) = 2*P1(:, 2:end-1);
            %phases
            P2 = angle(Y/L);
            P1 = P2(:, 1:floor(L/2)+1);
            f = Fs*(0:(L/2))/L;
            %all
            P2 = Y/L;
            P1 = P2(:, 1:floor(L/2)+1);
            f = Fs*(0:(L/2))/L;

        end

        function plot_time(data, t)
            if ~exist('max_value', 'var')
                max_value = max(data);
            end
            if ~exist('min_value', 'var')
                min_value = min(data);
            end
            range = max_value-min_value;

            plot(t, data , 'LineWidth', 2);
            xlabel('t (ms)')
            ylabel('PD')
%             xlim([0,1800])
            ylim([min_value-range*0.1, max_value+range*0.1])

        end
        function plot_frec(data, f, max_value)
            if ~exist('max_value', 'var')
                max_value = max(data)*1.1;
            end
            bar(f, data, 'LineWidth', 1) 

            ylim([0, max_value])

            xlabel('f (Hz)')
            ylabel('|P(f)|')
        end

    end
end
