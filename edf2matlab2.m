    function output = edf2matlab2(full_edf_name, output_folder_name, log, events2, vars2)
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
    clc
    tic;
    [file_path, file_name, ext] = fileparts(full_edf_name);
    if ~strcmpi(ext, '.edf')
        return;
    end
    mat_path = 'mat_files';
    if ~exist(mat_path, 'dir')
        mkdir(mat_path);
    end
    print_log(['Initiate loading and converting EDF file: ' strrep(file_name, '_', '\_') ext], log);
    full_edf_name_fixed = strrep(full_edf_name, ' ', '\ ');
    
    if length(full_edf_name)>100
        print_log('Error: the length of the EDF path exceeds the limit of 100 characters', log);
        return;
    end
    
    if isunix && ~ismac
        mat_file_path = strcat(file_path, filesep, mat_path, filesep, file_name, '_extend.mat');
        if(~exist(mat_file_path, 'file'))
            full_path = mfilename('fullpath');
            
            path = strrep(full_path(1:end-11), ' ', '\ ');
            setenv('LD_LIBRARY_PATH', [matlabroot '/bin/glnxa64:' matlabroot '/sys/os/glnxa64'])
            convertor = 'edf2mat/edfmat_ubuntu32';
            if ~isempty(strfind(computer, '64'))
                convertor = 'edf2mat/edfmat_ubuntu64';
            end
            system([path convertor, ' "', full_edf_name, '" "', mat_file_path, '"'])
            print_log(['Loading mat file: ' strrep(file_name, '_', '\_') '.mat'], log);
            clc
        end
        input_data = load(mat_file_path);

        input_data = input_data.Edf2Mat;
%         delete([output_folder_name filesep file_name '_extend.mat']);
    end
    

    if ispc
        input_data = edfmex(full_edf_name);
    end
    if ismac
        input_data_mac = Edf2Mat(full_edf_name_fixed);
        input_data = input_data_mac.RawEdf;
    end
    clc
    
    print_log(['Loading data has been completed: ' num2str(toc) ' seconds'], log);    
    
    %% Get all the indexes of user's messages (scan events's table and note all the ids of user's messages - 28 / 50)

    events      = input_data.FEVENT;
    parsedby_id = 0;
    tic;
    print_log('Initiate messages loading', log);    

    for i= 1: size(events,2)
        if strcmp(events(i).codestring, 'MESSAGEEVENT')
            parsedby_id = parsedby_id + 1;
            parsedby(1, parsedby_id) = i;
        end
    end

    
    %% Get all the user's messages (and the relevant timestems)
    tic;
    event_msgs      = cell(parsedby_id, 1);
    event_timestamps = zeros(1, parsedby_id);
    parsedby_id = 0;
    for i = parsedby
        parsedby_id = parsedby_id + 1;
        event_timestamps(parsedby_id) = events(i).sttime;
        event_msgs{parsedby_id}       = events(i).message;
    end
    
    print_log(['Loading messages has been completed: ' num2str(toc) ' seconds'], log);    

    %%  Get gase's information & timestamps
    tic;
    print_log('Initiate pupil data loading', log);    

    samples     = input_data.FSAMPLE;
    pupil_sizes = samples.pa;
    pupil_x     = samples.px;
    pupil_y     = samples.py;
    timestamps  = samples.time;
    pupil_size  = pupil_sizes(1,:);

    if pupil_size(1)<0  % Take "good" pupil (value >=0) 
        pupil_size = pupil_sizes(2,:);
    end
    if(pupil_sizes(1)>=0 && pupil_sizes(2)>=0)
        pupil_sizes(pupil_sizes==0) = nan;
        pupil_size = mean(pupil_sizes, 'omitnan');
    end
%     pupil_size = -diff(pupil_x);
    data.pupil_x    = pupil_x';
    data.pupil_y    = pupil_y';
    data.pupil_size = pupil_size';
    data.timestamps = timestamps';
    print_log(['Loading pupil data has been completed ' num2str(toc) ' seconds'], log);
    %% Get variables, messages & trial ids
    tic;
    
    %% find rate
    mode_id     = find(~cellfun(@isempty, strfind(event_msgs, '!MODE')), 1);
    mode_data   = strsplit(event_msgs{mode_id});
    data.rate   = str2double(char((mode_data(4)))); 
    data.file_name  = full_edf_name; 
    
    timestamps = timestamps/(1000/data.rate);
    event_timestamps = event_timestamps/(1000/data.rate);

    data.timestamps = data.timestamps/(1000/data.rate);
  
    %% find trials
    print_log('Parsing trials', log); 
    
    trial_ids = find(~cellfun(@isempty, strfind(event_msgs, 'TRIALID')));

  
    if(isempty(trial_ids))
        print_log('Error: Trials were not found', log);
        return;
    end
    trial_data.trial_names     = cellfun(@(x) str2double(char(regexp(char(x),'\d+','match'))), event_msgs(trial_ids));
    
    trial_data.Trial_Onset_num = arrayfun(@(x) get_Trial_Onset(x, event_timestamps, timestamps), trial_ids);
    
    Trial_Offset_ids = find(~cellfun(@isempty, strfind(event_msgs,'TRIAL_END')));

    if(isempty(Trial_Offset_ids))
        Trial_Offset_ids = find(~cellfun(@isempty, strfind(event_msgs,'TRIAL_RESULT 0')));
    end

    if(isempty(Trial_Offset_ids))
        Trial_Offset_ids              = vertcat(trial_ids(2:end), -1);
        trial_data.Trial_Offset_num   = arrayfun(@(x) get_Trial_Offset(x, event_timestamps, timestamps), Trial_Offset_ids)-1;
    else
        trial_data.Trial_Offset_num = arrayfun(@(x) get_Trial_Offset_msg(x, event_timestamps, timestamps), Trial_Offset_ids);
    end
    
    trial_data.trial_length    = trial_data.Trial_Offset_num-trial_data.Trial_Onset_num;
        
    data.trial_data = struct2table(trial_data);
    data.total_var_data_table = [];
    %% find variables
    print_log('Parsing variables', log);    
   
    data.total_var_data_table = [];
    var_ids  = find(~cellfun(@isempty, strfind(event_msgs,'!V TRIAL_VAR ')));

    if ~isempty(var_ids)
        total_var_data = parse_data.parse_vars(event_msgs, event_timestamps, event_timestamps, var_ids, trial_ids);
        data.total_var_data_table = struct2table(total_var_data);
    end
    
    external_var_ids = [];
    var_list_file = strcat(file_path, filesep, 'var_list.csv');
    if exist(var_list_file, 'file') 
        
        if exist('detectImportOptions', 'file')
            var_file = readtable(var_list_file, detectImportOptions(var_list_file));
        else
            var_file = readtable(var_list_file, 'Delimiter', ',');
        end
        
        vars2    = var_file.var_name;
        for i=1:size(vars2, 1)
            expression = strcat('\w*', char(vars2(i)), '*');
            external_var_ids = find(~cellfun(@isempty, regexp(event_msgs, expression, 'ONCE'))); 
        end        
    end
    
    if ~isempty(external_var_ids)
        total_var_data = parse_data.parse_vars(event_msgs, event_timestamps, event_timestamps, external_var_ids, trial_ids, true);
        data.total_var_data_table = struct2table(total_var_data);
    end

    vars_file = strcat(file_path, filesep, file_name, '_vars.csv');
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
    
    event_names = [];
    event_ids = [];
%     alternative_events = find(~cellfun(@isempty, regexp(event_msgs, '^\d *.','match')));
%     for event = 1:length(alternative_events)
%         event_arr = event_msgs(alternative_events(event));
%         event_data_arr  = strsplit(cell2mat(event_arr));
%         
%         if length(event_data_arr)==2
%             event_ids = [event_ids; alternative_events(event)];
%             if ~strcmp(event_names, event_data_arr{2});
%                 event_names = [event_names; event_data_arr(2)];
%             end
%         end
%     end
    if(isempty(event_ids))
        event_ids = find(~cellfun(@isempty, strfind(event_msgs,'!E TRIAL_EVENT_VAR')));
    end
    if ~isempty(event_ids)
        event_full_data = event_msgs(event_ids);
        event_full_timestamps = event_timestamps(event_ids);

        event_data_table = parse_data.parse_events(event_full_data, event_full_timestamps, timestamps, data.trial_data.Trial_Onset_num, data.trial_data.Trial_Offset_num, ~isempty(event_names));
        data.total_var_data_table = [data.total_var_data_table, event_data_table];
    else
        if(isempty(events2)) % events with general format - same to all participants
            [csv_file_name, csv_path_name]  = uigetfile('*.csv','Select the CSV file for events');
            full_csv_name                   = [csv_path_name csv_file_name];
            if(csv_path_name)
                data.full_csv_name = full_csv_name; 
                
                event_file = readtable(full_csv_name, 'Delimiter', ',');
                try
                    events2 = event_file.event_name;
                catch
                    print_log('Error: Wrong event file! ', log);
                    return;
                end
            end
        end
        if(~isempty(events2))
            for i=1:length(events2)
                event_ids_one = strcmp(event_msgs, events2{i}); 
                if(sum(event_ids_one)==0)
                    event_ids_one =~cellfun(@isempty, regexp(event_msgs, ['\d+ ', events2{i}, '$'], 'match')); 
                end
                event_msgs(event_ids_one, :) = {['!E TRIAL_EVENT_VAR ' strrep(strrep(events2{i}, ':', ''), ' ', '_')]};
            end            
        end
        event_ids = find(~cellfun(@isempty, strfind(event_msgs,'!E TRIAL_EVENT_VAR')));
        if ~isempty(event_ids)
            event_full_data = event_msgs(event_ids);
            event_full_timestamps = event_timestamps(event_ids);
            event_data_table = parse_data.parse_events(event_full_data, event_full_timestamps, timestamps, data.trial_data.Trial_Onset_num, data.trial_data.Trial_Offset_num, ~isempty(event_names));
            data.total_var_data_table = [data.total_var_data_table, event_data_table];
        end
    end
    
    data.events2 = events2;
    data.vars2   = vars2;
    trial_data.trial_length = trial_data.Trial_Offset_num-trial_data.Trial_Onset_num;
    if length(trial_data.trial_length) ~= size(data.total_var_data_table, 1)
        print_log('Error: The number of variables does not match the number of trials', log);
        return;
    end
    data.total_var_data_table.event_Trial_Offset = trial_data.trial_length;
    
    mm_file = strcat(file_path, filesep, file_name, '_mm.csv');
    if exist(mm_file, 'file')
        if exist('detectImportOptions', 'file')
            ap_data = readtable(mm_file, detectImportOptions(mm_file));
        else
            ap_data = readtable(mm_file, 'Delimiter', ',');
        end

        ap_mm     = ap_data.mm;
        ap_pixels = ap_data.pixels;
        
        pupil_data_type_id  = find(~cellfun(@isempty, strfind(event_msgs, 'PUPIL_DATA_TYPE')), 1);
        pupil_data_type_arr = strsplit(event_msgs{pupil_data_type_id});

        if (strcmp(pupil_data_type_arr(2), 'RAW_AUTOSLIP')) % measure by area
            ratio                 = ap_mm/(2*sqrt(ap_pixels/pi));
            pupil_diameter_pixels = 2*sqrt(data.pupil_size/pi);
        else
            ratio = ap_mm/ap_pixels;
            pupil_diameter_pixels = data.pupil_size;
        end
        
        pupil_diameter_mm     = pupil_diameter_pixels*ratio;
        data.pupil_size       = pupil_diameter_mm;
    end
    
    save([output_folder_name filesep file_name '.chp'], 'data');
    output = data;
    return;
end



function Trial_Onset_id = get_Trial_Onset(index, event_timestamps, timestamps)
    start_timestamp = event_timestamps(index);
    Trial_Onset_id  = find(timestamps>=start_timestamp, 1, 'first');
    if(isempty(Trial_Onset_id))
        Trial_Onset_id = 1;
    end
end


function Trial_offset_id = get_Trial_Offset_msg(index, event_timestamps, timestamps)
    start_timestamp = event_timestamps(index);
    Trial_offset_id  = find(timestamps>=start_timestamp, 1, 'first');
    if(isempty(Trial_offset_id))
        Trial_offset_id = length(timestamps);
    end
end


function Trial_Offset_id = get_Trial_Offset(index, event_timestamps, timestamps)
    next_timestamp = double(timestamps(end)+1000);
    if(index>0)
        next_timestamp = event_timestamps(index);
    end
    Trial_Offset_id = find(timestamps<next_timestamp,1,'last');
end

