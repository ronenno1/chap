function [pupil_data, blinks_data_positions] = fix_blinks2(pupil_data, Zoutliers, zerosOutliers, rate, linear_interpulation, gradient, debug_mode)
    pupil_data2 = pupil_data; 
    blinks_data_positions = [];
    gap_interval          = 100;                             % set the interval between two sets that appear consecutively for concatenation.
    gap_interval_samples  = gap_interval/(1000/rate); 
    start_debug = 0;
    if (~exist('debug_mode', 'var'))
        debug_mode = false;
    else
        debug_mode;
    end
    if(debug_mode && debug_mode>start_debug)
        pd =  pupil_data;
        pd(pd==0)=nan;    

        fig = figure(10);

        hold off
        plot(pd, 'black', 'LineWidth',2)
        hold on
    end
    
    pupil_data(pupil_data==0)=nan;    
    pupil_mean = mean(pupil_data, 'omitnan');
    pupil_std = std(pupil_data, 'omitnan');
    opd = pupil_data;
    if(Zoutliers>0)
        maxpupil2remove = Zoutliers*pupil_std+pupil_mean;
        pupil_data(pupil_data>=maxpupil2remove)=nan;
        minpupil2remove = -Zoutliers*pupil_std+pupil_mean;
        pupil_data(pupil_data<=minpupil2remove)=nan;
    end
    pupil_data(isnan(pupil_data))=0;    

    if(debug_mode && debug_mode>start_debug)
        pd =  pupil_data;
        pd(pd==0)=nan;    
        fig = figure(10);

        hold off
        plot(pd,'red', 'LineWidth',3)
        title(num2str(debug_mode))
        plot(opd,'yellow', 'LineWidth',3)

        hold on
    end
    num_of_zeros = sum(pupil_data==0);
    if(zerosOutliers>0 && num_of_zeros/size(pupil_data,1)*100>=zerosOutliers)
        pupil_data(pupil_data>=0)=0;
    end
    samples2smooth = ceil(rate/100);
    try
        smooth_data = smooth(pupil_data, samples2smooth, 'moving');
    catch
        smooth_data = my_smooth(pupil_data, samples2smooth, 'moving');
    end
    if(debug_mode && debug_mode>start_debug)
        pd_smooth =  smooth_data;
        pd_smooth(pd_smooth==0)=nan;    
        plot(pd_smooth,'black', 'LineWidth',3)
        hold on
        
%         t=1:size(smooth_data,1);
%         Fs = rate; % Sampling frequency
%         
%         Y = fft(pupil_data); % FFT of time domain signal y(t)
%         L = length(Y); % Length of Y
%         f = Fs*linspace(0,1,L); % Frequency vector (Hz)
% 
%         cutoff = 100; % cut off frequency (Hz)
%         NY = Y;
%         NY(f>cutoff & f<(Fs-cutoff)) = 0; 
%         Ny = real(ifft(NY));
%         hold off
% 
%         plot(Ny, 'r')
%         hold on
%         pk = pupil_data;
%         pk(pk==0)=nan;
%         plot(pk, 'b')

        
    end
    
    % pupil_data==0 return matrix of zeros and ones one the pupil is 0
    
    % find(diff(pupil_data==0)==1) return the first samples before the
    % missing values: diff = p(n+1)-p(n)
    
    % find(diff(pupil_data==0)==-1) return the last samples after the
    % missing values : diff = p(n+1)-p(n)
    
    % blinks look like: 0000001111110000
    % and diff : 0000001000000-1000 
    blinks      = vertcat(-1.*find(diff(pupil_data==0)==1), find(diff(pupil_data==0)==-1)+1);    
   
    without_noise = 0;
    gradient_approach = 0;
    diff_data = diff(pupil_data);
    pupil_data2(pupil_data2==0)=nan;
    diff_data2 = diff(pupil_data2);

    diff_mean = mean(diff_data, 'omitnan');
    diff_std  = std(diff_data, 'omitnan');
    
%     gradient    = diff_mean+gradient_crit*diff_std;
    blinks_data = -1.*(diff_data <=-gradient) + (diff_data >=gradient);
    if(size(blinks, 1)>0 && blinks(end)<0) && pupil_data(end)==0 
        blinks = vertcat(blinks, size(pupil_data, 1));
    end

    if gradient_approach>0
        blinks      = vertcat(-1.*find(blinks_data == -1), find(blinks_data == 1));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    [~,idx]     = sort(abs(blinks));
    blinks      = blinks(idx);

    % the trail starts with blink
    if(size(blinks, 1)>0 && blinks(1)>0) && pupil_data(1)==0 
        blinks = vertcat(-1, blinks);
    end
    
    % the trail ends with blink
    if(size(blinks, 1)>0 && blinks(end)<0) && pupil_data(end)==0 
        blinks = vertcat(blinks, size(pupil_data, 1));
    end

    % we doesn't have blinks!
    if(size(blinks)<2)          
        return;
    end
    smooth_data(smooth_data==0)=nan;
    diff_smooth_data    = diff(smooth_data);
    blink = 1;
    prev_p2 = -1;
    prev_p3 = -1;
    while blink < size(blinks, 1)
        while blink  <= size(blinks,1) && blinks(blink)>0  
            blink = blink + 1;
        end
        
        if blink > size(blinks)
            break;
        end

        blink_start = blinks(blink);
        
        blink = blink + 1; 
        while blink  <= size(blinks, 1) && blinks(blink)<0  
            blink = blink + 1;
        end
        while blink  < size(blinks, 1) && blinks(blink+1)>=0 && blinks(blink+1) == blinks(blink)+1   
            blink = blink + 1;
        end
        if blink > size(blinks)
            pupil_data(abs(blinks(blink-1)):end) = pupil_data(abs(blinks(blink-1)));
            if(debug_mode && debug_mode>start_debug)
                p2 = abs(blink_start);
                plot(p2, pd(p2), 'r*', 'LineWidth', 6);
            end
            break;
        end
        
        blink_end = blinks(blink);
        while blink  < size(blinks,1) && blinks(blink+1)>=0  
            blink = blink + 1;
        end
        strcat(num2str(blink_start),' :  ', num2str(blink_end-1));

        short2  = diff_smooth_data(2:abs(blink_start));
        
        short3  = diff_smooth_data(blink_end:end);
        
        p2      = find(short2>0, 1,'last'); %abs(blink_start);
          
        if isempty(p2)
            p2 = 2;
            blinks_data_positions = [blinks_data_positions, 0];
        else
            blinks_data_positions = [blinks_data_positions, -p2];
        end

        if pupil_data(p2+2)>0
            p2 = p2+2;
            blinks_data_positions(end) = -p2;

        end

        
        if(prev_p3>0 && p2-gap_interval_samples<=prev_p3)
            p2 = prev_p2;
            blinks_data_positions(end) = -prev_p3;

        end
        
        p3      = blink_end+find(short3<0, 1);
        if(isempty(p3))
            pupil_data(p2+1:end) = pupil_data(p2);            
            if(debug_mode && debug_mode>start_debug)
                plot(p2, pd(p2), 'r*', 'LineWidth',6);
            end
            p3 = size(pupil_data, 1)+1;
        end
        
        if pupil_data(p3-1)>0 || p3 == size(pupil_data, 1)+1
             p3 = p3-1;
        end
        
        blinks_data_positions = [blinks_data_positions, p3];
       
        if without_noise
            if(p2>2)
                p2 = abs(blink_start);
            end
            p3 = blink_end+1;
            [p2, p3];
        end
%         blinks_data_positions = [blinks_data_positions, [p2, p3]];
%         if without_noise
%             p2 = max(1, abs(blink_start)-100);
%             p3 = min(blink_end+1+100, size(pupil_data,1));
%             [p2,p3];
%         end
        p1      = max(1,p2-floor((p3-p2)/2));
        p4      = min(p3+floor((p3-p2)/2), size(pupil_data,1));
        if(p4==p3)
            p3 = p3-1;
        end
        if(debug_mode && debug_mode>start_debug)

            plot(p2, pd(p2), 'c*', 'LineWidth',10);
            plot(p3, pd(p3), 'm*', 'LineWidth',10);
            
        end

        if p2==2 || p2>p3
            pupil_data(1:p3) =  pupil_data(p3);
            continue;
        end
        if(p3>=p4)
            continue;
        end
        if(p2==p3) 
            continue;
        end
        if(linear_interpulation || p2-p1~=p4-p3)
            b = spline([p2 p3],[pupil_data(p2) pupil_data(p3)], p2:p3)';
        else
            b = spline([p1 p2 p3 p4],[pupil_data(p1) pupil_data(p2) pupil_data(p3) pupil_data(p4)],p2:p3)';
        end
        pupil_data = vertcat(pupil_data(1:p2-1), b, pupil_data(p3+1:end));
        prev_p2 = p2;
        prev_p3 = p3;
    end
  
    if(debug_mode && debug_mode>start_debug && ~isempty(blinks))
        plot(pd,'red', 'LineWidth',3)
        legend('original','smooth', 'blink onset', 'blink offset')
     
        set(fig, 'color', 'white');
        set(gca,'FontWeight','bold');
        present_threshold = false;
        if(present_threshold)
            diff_data = diff_data2;
            plot(diff_data, 'LineWidth',2);
            
            plot(diff_smooth_data, 'g', 'LineWidth',2);

            set(fig, 'color', 'white');
            set(gca,'FontWeight','bold');
            xL = get(gca,'XLim');
            line(xL,[gradient gradient],'Color','r');
            line(xL,[0 0],'Color','r');
            line(xL,[-gradient -gradient],'Color','r');
        end
        plot(pupil_data, 'LineWidth',2);
        xt = get(gca, 'XTick');                            
        set(gca, 'XTick', xt, 'XTickLabel', xt*(1000/rate));
        xlabel({'Time (ms)'})

        key = 0;

        while ~key
            key = waitforbuttonpress;
        end
        close(10);
        blinks_data_positions = blinks_data_positions*(1000/rate);
    end
    
    
    id = 1;
    while(id<size(blinks_data_positions, 2)-1)
        if(blinks_data_positions(id)>0 && blinks_data_positions(id)==-blinks_data_positions(id+1))
            blinks_data_positions(id:id+1) = [];
        else
            id = id+1;
        end
    end
    blinks_data_positions = abs(blinks_data_positions);
%     [n, bin] = histc(blinks_data_positions, unique(blinks_data_positions));
%     multiple = find(n > 1);
%     index    = find(ismember(bin, multiple));
%     blinks_data_positions(index) = [];

end