function blinks_data_positions = based_noise_blinks_detection(pupil_data, sampling_rate_in_hz) 
    
    blinks_data_positions = [];
    sampling_interval     = 1000/sampling_rate_in_hz; % compute the sampling time interval in milliseconds.
    gap_interval          = 100;                      % set the interval between two sets that appear consecutively for concatenation.
    
    %% Setting the blinks' candidates array
    % explanations for line 17:
    % pupil_data==0 returns a matrix of zeros and ones, where one means missing values for the pupil (missing values represented by zeros).
    % it looks like: 0000001111110000
    % diff(n) = pupil_data(n+1)-pupil_data(n)
    % find(diff(pupil_data==0)==1) returns the first sample before the missing values 
    % find(diff(pupil_data==0)==-1) returns the last missing values
    % it looks like: 00000100000-1000 
    % blink onset is represented by a negative value and blink offset is represented by a positive value
    blinks      = vertcat(-1.*find(diff(pupil_data==0)==1), find(diff(pupil_data==0)==-1)+1);    
    
    % Case 1: there are no blinks
    if(isempty(blinks))          
        return;
    end;
    
    % Sort the blinks by absolute value. in this way we are getting an array of blinks when the offset appears after the onset 
    [~, idx] = sort(abs(blinks));
    blinks   = blinks(idx);

    %% Edge cases
    % Case 2: the data starts with a blink. In this case, blink onset will be defined as the first missing value.
    if(size(blinks, 1)>0 && blinks(1)>0) && pupil_data(1)==0 
        blinks = vertcat(-1, blinks);
    end;
    
    % Case 3: the data ends with a blink. In this case, blink offset will be defined as the last missing sample
    if(size(blinks, 1)>0 && blinks(end)<0) && pupil_data(end)==0 
        blinks = vertcat(blinks, size(pupil_data, 1));
    end;

    %% Smoothing the data in order to increase the difference between the measurement noise and the eyelid signal.
    ms_4_smooting  = 10;                                    % using a gap of 10 ms for the smoothing
    samples2smooth = ceil(ms_4_smooting/sampling_interval); % amount of samples to smooth 
    smooth_data    = smooth(pupil_data, samples2smooth);    

    smooth_data(smooth_data==0) = nan;                      % replace zeros with NaN values
    diff_smooth_data            = diff(smooth_data);
    
    %% Finding the blinks' onset and offset
    blink                 = 1;                         % initialize blink index for iteration
    blinks_data_positions = zeros(size(blinks, 1), 1); % initialize the array of blinks
    prev_offset           = -1;                        % initialize the previous blink offset (in order to detect consecutive sets)    
    while blink < size(blinks, 1)
        % set the onset candidate
        onset_candidate = blinks(blink);
        blink           = blink + 1;  % increase the value for the offset
        
        % set the offset candidate
        offset_candidate = blinks(blink);
        blink            = blink + 1;  % increase the value for the next blink
        
        % find blink onset
        data_before = diff_smooth_data(2:abs(onset_candidate)); % returns all the data before the candidate
        blink_onset = find(data_before>0, 1, 'last');           % returns the last 2 samples before the decline
        
        % Case 2 (the data starts with a blink. In this case, blink onset will be defined as the first missing value.)
        if isempty(blink_onset)
            blink_onset = abs(onset_candidate);
        end;
        
        % correct the onset if we are not in case 2
        if blink_onset>0 || pupil_data(blink_onset+2)>0
            blink_onset      = blink_onset+2;
        end;
        
        % find blink offset
        data_after  = diff_smooth_data(abs(offset_candidate):end); % returns all data after the candidate
        blink_offset = offset_candidate+find(data_after<0, 1);     % returns the last sample before the pupil increase

        % Case 3 (the data ends with a blink. In this case, blink offset will be defined as the last missing sample.)
        if isempty(blink_offset)
            blink_offset = size(pupil_data, 1);
        end;
        
        % Set the onset to be equal to the previous offset in case where several sets of missing values are presented consecutively
        if (sampling_interval*blink_onset > gap_interval && sampling_interval*blink_onset-sampling_interval*prev_offset<=gap_interval)
            blink_onset = prev_offset;
        end;
        
        prev_offset = blink_offset-1;
        % insert the onset into the result array
        blinks_data_positions(blink-2) = sampling_interval*blink_onset;
        % insert the offset into the result array
        blinks_data_positions(blink-1) = sampling_interval*(blink_offset-1);
    end;
    
    %% Removing duplications (in case of consecutive sets): [a, b, b, c] => [a, c]
    [n, bin] = histc(blinks_data_positions, unique(blinks_data_positions));
    multiple = find(n > 1);
    index    = find(ismember(bin, multiple));
    blinks_data_positions(index) = [];
end