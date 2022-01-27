function merge_chps(input_folder_name, subject_id)
    if ~exist('input_folder_name', 'var')
        input_folder_name = uigetdir('select directory');
        if(~input_folder_name)
            return;
        end
    end
    if ~exist('subject_id', 'var')
        subject_id  = ''; 
        prompt1         = {'Subject id'};
        dlg_title1      = 'Subject id ';
        num_lines       = 1;
        answer1         = inputdlg(prompt1, dlg_title1, num_lines, {subject_id});
        if sum(size(answer1))>0 && ~strcmp(answer1{1}, '')
            subject_id  = answer1{1};
        else
            disp('Error: Condition names have to be chosen');
            return;
        end
    end
    file_names = dir([input_folder_name filesep '*.chp']);
    file_names = {file_names.name}';
    
    for file_id = 1:length(file_names)
        files{file_id} = load([input_folder_name filesep file_names{file_id}], '-mat');
        times(file_id) = files{file_id}.data.timestamps(1,:);
    end
    [~, sort_ids] = sort(times);
    files = files(sort_ids);
    output_data = files{1};
    [~, ~, ext] = fileparts(output_data.data.file_name);

    output_data.data.file_name = [subject_id, ext];

    for file_id = 2:length(file_names)
        data_file = files{file_id}.data;
        number_of_samples = length(output_data.data.pupil_size);

        data_file.trial_data.Trial_Onset_num  = data_file.trial_data.Trial_Onset_num+number_of_samples;
        data_file.trial_data.Trial_Offset_num = data_file.trial_data.Trial_Offset_num+number_of_samples;
        data_file.trial_data.trial_names      = data_file.trial_data.trial_names+length(output_data.data.trial_data.trial_names);
        field_names = fieldnames(data_file);
        for field = 1:length(field_names)
            if strcmp(field_names{field}, 'rate') || strcmp(field_names{field}, 'file_name') || strcmp(field_names{field}, 'full_csv_name')
                continue;
            end
            try
                output_data.data.(field_names{field}) = [output_data.data.(field_names{field}); data_file.(field_names{field})];
            catch
                disp('Error: All files have to include the same variables and events!');
                return;
            end
        end
    end

    data = output_data.data;
    save([input_folder_name, filesep, subject_id, '.chp'], 'data');
    
    data_file =  [input_folder_name, filesep, subject_id, ext];
      
    fileID = fopen(data_file,'w');
    fprintf(fileID, '%s ', 1:5000);
    fclose(fileID);
end