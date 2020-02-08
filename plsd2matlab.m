function output = plsd2matlab(full_file_name, output_folder_name, log, events2, vars2)
    
    if (~exist('log', 'var'))
        log = '';
    end
    if (~exist('events2', 'var'))
        events2 = [];
    end    
    if (~exist('vars2', 'var'))
        vars2 = [];
    end
    output = [];
    [path, file_name, ext] = fileparts(full_file_name);
    if ~strcmpi(ext, '.plsd')
        return;
    end
      
    
    info_csv_name = [path filesep file_name '_info.csv'];
    if ~exist(info_csv_name, 'file')
        print_log(strcat('Error (3): info file does not found, please add: ', strrep(strcat(file_name, '_events.csv') ,'_','\_') , ' to ', path), log);    
        return;
    end

    event_csv_name = [path filesep file_name '_events.csv'];
    if ~exist(event_csv_name, 'file')
        print_log(strcat('Error (3): events file does not found, please add: ', strrep(strcat(file_name, '_events.csv') ,'_','\_') , ' to ', path), log);    
        return;
    end


    
    print_log(['Start load and convert PLSD file: ' strrep(file_name, '_', '\_') ext], log);

    loaded = false;
    try
        
        if exist('detectImportOptions', 'file')
            info_file = readtable(info_csv_name, detectImportOptions(info_csv_name));
        else
            info_file = readtable(info_csv_name, 'Delimiter', ',');
        end     

        
        diff = str2num(char(info_file{5,2}))-str2num(char(info_file{6,2}));
        dat_file_name = [path filesep file_name '.dat'];

        copyfile(full_file_name, dat_file_name);
            
        
        if exist('detectImportOptions', 'file')
            pupil_data = readtable(dat_file_name, detectImportOptions(dat_file_name));
        else
            pupil_data = readtable(dat_file_name, 'Delimiter', ',');
        end

        % CHAP removes all the irrelevant fields (CHAP needs only 3)
        if (size(pupil_data, 2)>3)
        
            pupil_data_short.world_timestamp = pupil_data.world_timestamp;
            pupil_data_short.eye_id          = pupil_data.eye_id;
            pupil_data_short.diameter_3d     = pupil_data.diameter_3d;
            new_pupil_data_table = struct2table(pupil_data_short);
            
            writetable(new_pupil_data_table, [path filesep file_name '.dat']);    
            copyfile([path filesep file_name '.dat'], [path filesep file_name '.plsd']);
        end
        delete([path filesep file_name '.dat']);

        pupil_data.world_timestamp = pupil_data.world_timestamp+diff;

        first_sample_timestamp = pupil_data.world_timestamp(1)/86400 + datenum('01-Jan-1970');
        record_start_time = datenum([char(info_file{3, 2}) ' ' char(info_file{4, 2})], 'dd.mm.yyyy HH:MM:SS');

        seconds_diff = (record_start_time - first_sample_timestamp)*86400;
        hours_diff = round(seconds_diff/3600);
        addition2times = hours_diff*3600;

        
        pupil_data_r = pupil_data(pupil_data.eye_id==0, :);
        pupil_data_l = pupil_data(pupil_data.eye_id==1, :);

        rate = 0;

        for sample=2:100
            rate = max(rate, round(1/(pupil_data_l.world_timestamp(sample)-pupil_data_l.world_timestamp(sample-1))));
        end
        rate = round(rate/(10))*10;
        rate_s = 1/rate;

        gaps_r = calculate_gaps(pupil_data_r.world_timestamp, rate_s);
        gaps_l = calculate_gaps(pupil_data_l.world_timestamp, rate_s);

        pupil_r = insert_gaps(pupil_data_r, gaps_r, rate_s);
        pupil_l = insert_gaps(pupil_data_l, gaps_l, rate_s);

        first_eye  = pupil_r;
        second_eye = pupil_l;
        if(pupil_r.world_timestamp(1)>=pupil_l.world_timestamp(1))
            first_eye  = pupil_l;
            second_eye = pupil_r;
        end

        equal_times = false;
        while ~equal_times
            if second_eye.world_timestamp(1)-first_eye.world_timestamp(1) < rate_s/2
                equal_times = true;
            end
            first_eye(1,:) = [];
        end
        num_of_samples = min(size(first_eye, 1), size(second_eye, 1))-1;
        first_eye      = first_eye(1:num_of_samples, :);
        second_eye     = second_eye(1:num_of_samples, :);
        pupil_data_correct.world_timestamp = first_eye.world_timestamp;

        first_eye.diameter_3d(isnan(first_eye.diameter_3d))   = second_eye.diameter_3d(isnan(first_eye.diameter_3d));
        second_eye.diameter_3d(isnan(second_eye.diameter_3d)) = first_eye.diameter_3d(isnan(second_eye.diameter_3d));

        pupil_data_correct.pupil_size = (first_eye.diameter_3d+second_eye.diameter_3d)/2;
        
        data.timestamps = pupil_data_correct.world_timestamp;
        
        data.pupil_size = pupil_data_correct.pupil_size;
        data.pupil_x = pupil_data_correct.pupil_size;
        data.pupil_y = pupil_data_correct.pupil_size;
        loaded = true;
    catch
        
    end        
    if ~loaded
        print_log('Error: incompetible file (2), please check your file', log);    
        return;
    end
    print_log('Start loading pupil data', log);  
    timestamps = rate*data.timestamps;
    
    data.rate = rate;
    data.file_name  = full_file_name;     
    
    %% Get all the indexes of user's messages
    
    tic;
    try
        if exist('detectImportOptions', 'file')
            mes_data_table = readtable(event_csv_name, detectImportOptions(event_csv_name));
        else
            mes_data_table = readtable(event_csv_name, 'Delimiter', ',');
        end     
    catch
        print_log('Error (4): incompetible file, please check your file', log);    
        return;
    end
    print_log('starting loading messages', log);    

    event_msgs        = mes_data_table.message;
    event_timestamps  = rate*(86400*(datenum(mes_data_table.timestamp, 'yyyy-mm-dd HH:MM:SS.FFF')- datenum('01-Jan-1970'))-addition2times);

    print_log(['Finished loading messages: ' num2str(toc) ' seconds'], log);    

    %%  Get gase's information & timestamps
    tic;
    
    %% find trials
    print_log('Parsing trials', log);
    
    trial_ids = find(~cellfun(@isempty, strfind(event_msgs,'TRIALID'))); %trial is defined by message with the form TRIALID [num_of_trial]
    if(isempty(trial_ids))
        print_log('Error: trials did not found', log);
        return;
    end
    trial_data.trial_names     = cellfun(@(x) str2double(char(regexp(char(x),'\d+','match'))), event_msgs(trial_ids));        
    trial_data.Trial_Onset_num = arrayfun(@(timestamp) get_trial_data_onset(timestamp, timestamps), event_timestamps(trial_ids));
    
    Trial_Offset_ids              = ~cellfun(@isempty, strfind(event_msgs, 'TRIAL_END')); %trial is defined by message with the form TRIALID [num_of_trial]
    trial_data.Trial_Offset_num   = arrayfun(@(x) get_trial_data_onset(x, timestamps), event_timestamps(Trial_Offset_ids));
   
    trial_data.trial_length    = trial_data.Trial_Offset_num-trial_data.Trial_Onset_num;
    data.trial_data = struct2table(trial_data);
    
    data.total_var_data_table = [];
    print_log('Parsing variables', log);    
    var_ids = find(~cellfun(@isempty, strfind(event_msgs,'!V TRIAL_VAR')));

    if ~isempty(var_ids)
        total_var_data = parse_data.parse_vars(event_msgs, event_timestamps, timestamps, var_ids, trial_data.Trial_Onset_num);
        data.total_var_data_table = struct2table(total_var_data);
    end

    vars_file = strcat(path, filesep, file_name, '_vars.csv');
    if exist(vars_file, 'file')
        try
            var_data_table = parse_data.parse_external_vars(vars_file);
            data.total_var_data_table = [data.total_var_data_table, var_data_table];
        catch err
            print_log(['Error: ' err.message], log);
            return;
        end
    end

    if isempty(data.total_var_data_table)
        print_log('There are no variables, dummy variable is created...', log);
        dummy.dummy = repmat({'dummy'},size(trial_ids));
        data.total_var_data_table = struct2table(dummy);
    end
    
    print_log('Parsing events', log);
    data.event_data = [];
    event_ids = find(~cellfun(@isempty, strfind(event_msgs,'!E TRIAL_EVENT_VAR')));
        event_full_data = event_msgs(event_ids);
        event_full_timestamps = event_timestamps(event_ids);

    if ~isempty(event_ids)
        event_data_table = parse_data.parse_events(event_full_data, event_full_timestamps, timestamps, data.trial_data.Trial_Onset_num);
        data.total_var_data_table = [data.total_var_data_table, event_data_table];
    end
    
    data.events2 = events2;
    data.vars2   = vars2;
    trial_data.trial_length = trial_data.Trial_Offset_num-trial_data.Trial_Onset_num;
    data.total_var_data_table.event_Trial_Offset = trial_data.trial_length;

    save([output_folder_name filesep file_name '.chp'], 'data');
    output = data;
end



function start_time = get_trial_data_onset(timestamp, timestamps)
    start_time      = find(timestamp>timestamps, 1, 'last');
end

function pupil_data_full = insert_gaps(pupil_data, num_of_gaps, rate_s)
    total_samples = num_of_gaps + size(pupil_data, 1);
    pupil_data_full.diameter_3d     = zeros(total_samples, 1);
    pupil_data_full.world_timestamp = zeros(total_samples, 1);

    sample = 1;
    for id =2:size(pupil_data, 1)
        pupil_data_full.diameter_3d(sample, :)     = pupil_data.diameter_3d(id-1);
        pupil_data_full.world_timestamp(sample, :) = pupil_data.world_timestamp(id-1);

        if pupil_data.world_timestamp(id)-pupil_data.world_timestamp(id-1) > rate_s
            num_of_samples = round((pupil_data.world_timestamp(id)-pupil_data.world_timestamp(id-1))/rate_s)-1;
            for inner_sample = 1:num_of_samples
                pupil_data_full.world_timestamp(sample + inner_sample, :) = pupil_data_full.world_timestamp(sample, :) + inner_sample*rate_s;
            end
            sample = sample + num_of_samples;
        end
        sample = sample + 1;
    end
    pupil_data_full.diameter_3d(pupil_data_full.diameter_3d == 0) = NaN;
    pupil_data_full = struct2table(pupil_data_full);
end


function gaps = calculate_gaps(world_timestamp, rate_s)
    gaps = 0;
    for sample = 2:size(world_timestamp, 1)
        if world_timestamp(sample)-world_timestamp(sample-1) > rate_s
            gaps = gaps + round((world_timestamp(sample)-world_timestamp(sample-1))/rate_s)-1;
        end
    end
end
