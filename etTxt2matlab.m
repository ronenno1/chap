function output = etTxt2matlab(full_txt_name, output_folder_name, log, events2, vars2)
    
    if (~exist('log', 'var'))
        log = false;
    end
    if (~exist('events2', 'var'))
        events2 = [];
    end    
    if (~exist('vars2', 'var'))
        vars2 = [];
    end
    output = [];
    [path, file_name, ext] = fileparts(full_txt_name);
    if ~strcmpi(ext, '.txt')
        return;
    end
      
    event_csv_name = [path filesep file_name '_events.csv'];
    if ~exist(event_csv_name, 'file')
        print_log(strcat('Error (3): The events file was not found. Please add: ', strrep(strcat(file_name, '_events.csv') ,'_','\_') , ' to ', path), log);    
        return;
    end

    print_log(['Initiate loading and converting TXT file: ' strrep(file_name, '_', '\_') ext], log);

    % if ubuntu doesn't supprt python (like in ubuntu 22.04+, fix a problem
    % with:
    % cd /full_path_to_matlab_mcrroot/v94/bin/glnxa64
    % mv libexpat.so.1 libexpat.so.1.NOFIND

    loaded = false;
    try
        csv_file_name = strcat(path, filesep, file_name, '.csv');
        json_file_name = strcat(path, filesep, file_name, '.json');
        if ~exist(csv_file_name, 'file')
            [folder, ~, ~] = fileparts(mfilename('fullpath'));
            if isunix
                    status = system(['python3 "' folder filesep 'json2csv.py" ', '"', strcat(path, filesep, file_name), '"']);
                    if status
                        system(['python "' folder filesep 'json2csv.py" ', '"', strcat(path, filesep, file_name), '"']);
                    end
                
            else
                system(['py "' folder filesep 'json2csv.py" ', '"', strcat(path, filesep, file_name), '"']);
            end
            if exist(json_file_name, 'file')
               delete(json_file_name);
            end
        end
        
        print_log(['Initiate loading and converting csv file: ' strrep(file_name, '_', '\_') '.csv'], log);
        
        if exist('detectImportOptions', 'file')
            raw_data_table = readtable(csv_file_name, detectImportOptions(csv_file_name));
        else
            raw_data_table = readtable(csv_file_name, 'Delimiter', ',');
        end

        data.timestamps = raw_data_table.values_frame_timestamp;
        raw_data_table.values_frame_lefteye_psize(raw_data_table.values_frame_lefteye_psize==0) = raw_data_table.values_frame_righteye_psize(raw_data_table.values_frame_lefteye_psize==0);
        raw_data_table.values_frame_righteye_psize(raw_data_table.values_frame_righteye_psize==0) = raw_data_table.values_frame_lefteye_psize(raw_data_table.values_frame_righteye_psize==0);
        
        data.pupil_size = (raw_data_table.values_frame_lefteye_psize + raw_data_table.values_frame_righteye_psize)/2;

        raw_data_table.values_frame_lefteye_avg_x(raw_data_table.values_frame_lefteye_avg_x==0) = raw_data_table.values_frame_righteye_avg_x(raw_data_table.values_frame_lefteye_avg_x==0);
        raw_data_table.values_frame_righteye_avg_x(raw_data_table.values_frame_righteye_avg_x==0) = raw_data_table.values_frame_lefteye_avg_x(raw_data_table.values_frame_righteye_avg_x==0);
        data.pupil_x    = (raw_data_table.values_frame_lefteye_avg_x + raw_data_table.values_frame_righteye_avg_x)/2;
        
        raw_data_table.values_frame_lefteye_avg_y(raw_data_table.values_frame_lefteye_avg_y==0) = raw_data_table.values_frame_righteye_avg_y(raw_data_table.values_frame_lefteye_avg_y==0);
        raw_data_table.values_frame_righteye_avg_y(raw_data_table.values_frame_righteye_avg_y==0) = raw_data_table.values_frame_lefteye_avg_y(raw_data_table.values_frame_righteye_avg_y==0);
        data.pupil_y    = (raw_data_table.values_frame_lefteye_avg_y + raw_data_table.values_frame_righteye_avg_y)/2;
%         data.pupil_size = abs(raw_data_table.values_frame_lefteye_avg_x - raw_data_table.values_frame_righteye_avg_x);
%         data.pupil_size = data.pupil_y; using 
        loaded = true;
    catch
    end        
    if ~loaded
        print_log('Error: Incompatible file (2), please check your file', log);    
        return;
    end
    print_log('Initiate pupil data loading', log);  
    try
        timestamps = 86400*(datenum(data.timestamps(:), 'yyyy-mm-dd HH:MM:SS.FFF')- datenum('01-Jan-1970'))';
    catch
        timestamps = 86400*(datenum(data.timestamps(:))- datenum('01-Jan-1970'))';
    end
    data.rate = 0;
    first_index = 2;
    while data.rate==0
        first_index = first_index+1;
        data.rate   = round(1000/(1000*(timestamps(first_index) - timestamps(first_index-1))), -1);
    end
    

    force_sampling_rate = 0;
    if force_sampling_rate>0 && data.rate>force_sampling_rate
        ratio = data.rate/force_sampling_rate;
        data.pupil_size = data.pupil_size(ratio:ratio:end, :);
        data.pupil_y    = data.pupil_y(ratio:ratio:end, :); 
        data.pupil_x    = data.pupil_x(ratio:ratio:end, :); 
        data.timestamps = data.timestamps(ratio:ratio:end, :); 
        timestamps      = timestamps(:, 1:ratio:end) ; 
        data.rate       = force_sampling_rate;
    end
    
    

    [timestamps, data] = parse_data.add_missing_samples(timestamps, data, log);
    
    data.file_name  = full_txt_name;     
    timestamps = timestamps*data.rate;
    
    %% Get all the indexes of user's messages
    
    tic;
    try
        if exist('detectImportOptions', 'file')
            mes_data_table = readtable(event_csv_name, detectImportOptions(event_csv_name));
        else
            mes_data_table = readtable(event_csv_name, 'Delimiter', ',');
        end     
    catch
        print_log('Error (4): Incompatible file, please check your file', log);    
        return;
    end
    print_log('Initiate messages loading', log);    

    event_msgs        = mes_data_table.message;
    try
        event_timestamps  = data.rate*86400*(datenum(mes_data_table.timestamp, 'yyyy-mm-dd HH:MM:SS.FFF')- datenum('01-Jan-1970'));
    catch
        event_timestamps  = data.rate*86400*(datenum(mes_data_table.timestamp)- datenum('01-Jan-1970'));
    end
    print_log(['Loading messages has been completed: ' num2str(toc) ' seconds'], log);    

    %%  Get gase's information & timestamps
    tic;
    
    %% find trials
    print_log('Parsing trials', log);
    
    trial_ids = find(~cellfun(@isempty, strfind(event_msgs,'TRIALID'))); %trial is defined by message with the form TRIALID [num_of_trial]
    if(isempty(trial_ids))
        print_log('Error: Trials were not found', log);
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
    mm_file = strcat(path, filesep, file_name, '_mm.csv');
    if exist(mm_file, 'file')
        if exist('detectImportOptions', 'file')
            ap_data = readtable(mm_file, detectImportOptions(mm_file));
        else
            ap_data = readtable(mm_file, 'Delimiter', ',');
        end
        ap_mm     = ap_data.mm;
        ap_pixels = ap_data.pixels;
        ratio     = ap_mm/(2*sqrt(ap_pixels/pi));
        
        pupil_diameter_pixels = 2*sqrt(data.pupil_size/pi);
        pupil_diameter_mm     = pupil_diameter_pixels*ratio;
        data.pupil_size       = pupil_diameter_mm;
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


