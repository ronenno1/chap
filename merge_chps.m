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

    for file_id = 1:length(file_names)-1
        number_of_samples = length(files{file_id}.data.pupil_size);
        [~, ~, ext]       = fileparts(files{file_id}.data.file_name);
        files{file_id}.data.file_name = [subject_id, ext];
        files{file_id+1}.data.trial_data.Trial_Onset_num  = files{file_id+1}.data.trial_data.Trial_Onset_num+number_of_samples;
        files{file_id+1}.data.trial_data.Trial_Offset_num = files{file_id+1}.data.trial_data.Trial_Offset_num+number_of_samples;
        files{file_id+1}.data.trial_data.trial_names      = files{file_id+1}.data.trial_data.trial_names+length(files{file_id}.data.trial_data.trial_names);
        field_names = fieldnames(files{file_id}.data);
        for field = 1:length(field_names)
            if strcmp(field_names{field}, 'rate') || strcmp(field_names{field}, 'file_name')
                files{file_id+1}.(field_names{field}) = files{file_id}.data.(field_names{field});
                continue;
            end
            files{file_id+1}.(field_names{field}) = [files{file_id}.data.(field_names{field}); files{file_id+1}.data.(field_names{field})];
        end
    end

    data = files{end};
    save([input_folder_name, filesep, subject_id, '.chp'], 'data');
    
    data_file =  [input_folder_name, filesep, subject_id, ext];
      
    fileID = fopen(data_file,'w');
    fprintf(fileID, '%s ', 1:5000);
    fclose(fileID);
end