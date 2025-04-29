function output = etTxt2csv(full_txt_name)
    clc
    
    output = [];
    [path, file_name, ext] = fileparts(full_txt_name);
    if ~strcmpi(ext, '.txt')
        return;
    end
      
    output_path = ['.' filesep 'mat_files'];
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end
    event_csv_name = [path filesep file_name '_events.csv'];
    if ~exist(event_csv_name, 'file')
        display(strcat('Error (3): events file does not found, please add ', strrep(file_name,'_','\_'), '_events.csv to ', path));    
        return;
    end
    
    loaded = false;
    try
        csv_file_name = strcat(path, filesep, file_name, '.csv');
        json_file_name = strcat(path, filesep, file_name, '.json');

        if ~exist(csv_file_name, 'file')
            origin_file_path = strrep(strcat(path, filesep, file_name), ' ', '\ ');
            [folder, ~, ~] = fileparts(mfilename('fullpath'));
            if isunix
                system(['python ' folder filesep 'json2csv.py ', origin_file_path]);
            else
                system(['py  ' folder  filesep 'json2csv.py ', origin_file_path]);
            end
            if exist(json_file_name, 'file')
               delete(json_file_name);
            end

        end
        
        display(['Start loading and convert csv file: ' strrep(file_name, '_', '\_') ext]);
        
    
        if exist('detectImportOptions', 'file')
            raw_data_table = readtable(csv_file_name, detectImportOptions(csv_file_name));
        else
            raw_data_table = readtable(csv_file_name, 'Delimiter', ',');
        end

        
        data.timestamps = raw_data_table.values_frame_timestamp;
        data.pupil_size = (raw_data_table.values_frame_lefteye_psize + raw_data_table.values_frame_righteye_psize)/2;
        data.pupil_x    = raw_data_table.values_frame_lefteye_avg_x;
        data.pupil_y    = raw_data_table.values_frame_lefteye_avg_y;
        raw_data_table2 = struct2table(data);

        loaded = true;
    catch
    end        
    display('Start loading pupil data');    
    timestamps = 86400*(datenum(data.timestamps(:))- datenum('01-Jan-1970'))';

    data.rate = 0;
    first_index = 1;
    while data.rate==0
        first_index = first_index+1;
        data.rate   = round(1000/(1000*(timestamps(first_index) - timestamps(first_index-1))), -1);
    end

    data.file_name  = full_txt_name;     

    %% Get all the indexes of user's messages
    
    tic;
    
    try
        if exist('detectImportOptions', 'file')
            mes_data_table = readtable(event_csv_name, detectImportOptions(event_csv_name));
        else
            mes_data_table = readtable(event_csv_name, 'Delimiter', ',');
        end

    catch
    
        display('Error (4): incompetible file, please check your file');    
        return;
    end
        
    event_msgs        = mes_data_table.message;
    event_timestamps  = 86400*(datenum(mes_data_table.timestamp)- datenum('01-Jan-1970'));
    event_dates = mes_data_table.timestamp;

    
    message = cell(size(data.pupil_size));
    mes = table(message);
    raw_data_table2 = [raw_data_table2 mes];
    events_table = {};
    for mes = 1:size(event_timestamps, 1)
        timestamp = event_timestamps(mes);
        event_msg = event_msgs(mes);
        event_date= event_dates(mes);
        pos      = find(timestamp>timestamps, 1, 'last');
        T1 = {event_date, 0, 0, 0, event_msg};
        events_table  = [events_table; T1];
    end
    
    raw_data_table2 = [raw_data_table2; cell2table(events_table, 'VariableNames',{'timestamps' 'pupil_size' 'pupil_x' 'pupil_y' 'message'})];
    
    writetable(raw_data_table2, [file_name '.dat']);

end

