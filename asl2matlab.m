function output = asl2matlab(full_asl_name, output_folder_name, log, events2, vars2)
    
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
    [path, file_name, ext] = fileparts(full_asl_name);
    if ~strcmpi(ext, '.asl')
        return;
    end
      
    event_csv_name = [path filesep file_name '_events.csv'];
    if ~exist(event_csv_name, 'file')
        print_log(strcat('Error (3): events file does not found, please add: ', strrep(strcat(file_name, '_events.csv') ,'_','\_') , ' to ', path), log);    
        return;
    end

    print_log(['Start load and convert ASL file: ' strrep(file_name, '_', '\_') ext], log);

        
    print_log(['Start loading and convert csv file: ' strrep(file_name, '_', '\_') ext], log);
    copyfile(full_asl_name,[path filesep file_name '.dat']);
    
    raw_data_table  = readtable([path filesep file_name '.dat'], 'Delimiter', ',');
    delete([path filesep file_name '.dat']);
    
    data.timestamps = raw_data_table.times;
    data.pupil_size = raw_data_table.pupil_size;
    data.pupil_x    = raw_data_table.pupil_x;
    data.pupil_y    = raw_data_table.pupil_y;

    print_log('Start loading pupil data', log);    

    data.rate = 60;

    data.file_name  = full_asl_name;     

    %% Get all the indexes of user's messages
    
    tic;
    
    try
        mes_data_table = readtable(event_csv_name, 'Delimiter', ',');    
    catch
    
        print_log('Error (4): incompetible file, please check your file', log);    
        return;
    end
        
    event_msgs        = mes_data_table.message;
    event_timestamps  = mes_data_table.times';

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
    trial_data.Trial_Onset_num = event_timestamps(trial_ids)';
    
    Trial_Offset_ids              = ~cellfun(@isempty, strfind(event_msgs, 'TRIAL_END')); %trial is defined by message with the form TRIALID [num_of_trial]
    trial_data.Trial_Offset_num   = event_timestamps(Trial_Offset_ids)';
   
    trial_data.trial_length    = trial_data.Trial_Offset_num-trial_data.Trial_Onset_num;
    
    data.trial_data = struct2table(trial_data);
    
    print_log('Parsing variables', log);    
    var_ids = find(~cellfun(@isempty, strfind(event_msgs,'!V TRIAL_VAR')));
    if ~isempty(var_ids)
        total_var_data = parse_data.parse_vars(event_msgs, event_timestamps, data.timestamps, var_ids, trial_data.Trial_Onset_num);
        data.total_var_data_table = struct2table(total_var_data);
    end

    print_log('Parsing events', log);    
    data.event_data = [];
    event_ids = find(~cellfun(@isempty, strfind(event_msgs,'!E TRIAL_EVENT_VAR')));
    event_full_data = event_msgs(event_ids);
    event_full_timestamps = event_timestamps(event_ids);

    if ~isempty(event_ids)
        event_data_table = parse_data.parse_events(event_full_data, event_full_timestamps, data.timestamps, data.trial_data.Trial_Onset_num);
        data.total_var_data_table = [data.total_var_data_table, event_data_table];
    end

    mm_file = strcat(path, filesep, file_name, '_mm.csv');
    if exist(mm_file, 'file')
        ap_data   = readtable(mm_file, 'Delimiter', ',');
        ap_mm     = ap_data.mm;
        ap_pixels = ap_data.pixels;
        ratio     = ap_mm/(2*sqrt(ap_pixels/pi));
        
        pupil_diameter_pixels = data.pupil_size;
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