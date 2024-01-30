function output = load_chap_data(compressed_data, log)
    if (~exist('log', 'var'))
        log = '';
    end
    tic
    data.pupil_x    = compressed_data.pupil_x;
    data.pupil_y    = compressed_data.pupil_y;
    data.pupil_size = compressed_data.pupil_size;
    data.timestamps = compressed_data.timestamps;
    timestamps      = data.timestamps;
  
    %% find trials
    percentages = zeros(1, 10);
    print_log('Initiate trials writing', log);    

    data_var.trial_id        = -1*ones(1, size(timestamps, 1));
    
    for i = 1:size(compressed_data.trial_data, 1)
        percentage = round(100*(i/size(compressed_data.trial_data, 1)));
        if percentage>0 && ~mod(percentage, 10) && percentages(percentage/10)==0 
            print_log(['Trials were written: ' num2str(percentage) '%'], log);    
            percentages(percentage/10) = 1;
        end
        data_var.trial_id(compressed_data.trial_data.Trial_Onset_num(i): compressed_data.trial_data.Trial_Offset_num(i)) = i;
    end
    
    %% find variables
    print_log('Initiate variables writing', log);    
    data_table     = compressed_data.total_var_data_table;
    percentages = zeros(1, 10);

    headers = data_table.Properties.VariableNames';
    var_names   = headers(cellfun(@isempty, strfind(headers,'event_')));
    event_names = headers(~cellfun(@isempty, strfind(headers,'event_')));
   
    for i = 1:size(var_names, 1)
        data_event2.(var_names{i}) = cell(1,size(timestamps, 1));
    end
   
    for i = 1:size(event_names, 1)
        data_event2.(event_names{i}) = zeros(1,size(timestamps, 1));
    end

    for trial = 1:size(compressed_data.trial_data, 1)
        percentage = round(100*(trial/size(compressed_data.trial_data, 1)));
        if percentage>0 && ~mod(percentage, 10) && percentages(percentage/10)==0 
            print_log(['Variables were written: ' num2str(percentage) '%'], log);    
            percentages(percentage/10) = 1;
        end

        for var = 1:size(headers, 1)
            var_name = headers(var);
            data_event2.(char(var_name))(data_var.trial_id == trial) = compressed_data.total_var_data_table.(char(var_name))(trial, :);
        end
    end

    var_namesk = fieldnames(data_var);
    for i = 1: size(var_namesk,1)
        v_name          = char(var_namesk(i));
        data.(v_name)   = data_var.(v_name)';
    end 

    
    for i = 1: size(headers, 1)
        e_name          = char(headers(i));
        data.(e_name)   = data_event2.(e_name)';
    end
    
    full_data_table       = struct2table(data);
    output.data           = full_data_table;

    output.var_data_table = compressed_data.total_var_data_table;
    output.rate           = compressed_data.rate;
    output.events2        = compressed_data.events2;
    output.vars2          = compressed_data.vars2;
    output.file_name      = compressed_data.file_name;
    
    print_log(['Table creation has been completed: ' num2str(toc) ' seconds'], log);    
    return;
end
