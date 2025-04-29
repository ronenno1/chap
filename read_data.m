classdef read_data
    methods(Static)
        
        function data = read_file(log) % read result file and convert it to chap format
            start = tic;
            data = [];
            [file_name, path_name] = uigetfile({'*.edf;*.txt;*.asl;*.tbi;*.plsd;*.dat',  'Eye-link(*.edf) / Eye-tribe(*.txt) / ASL(*.asl) / Tobii(*.tbi) / pupil-labs(*.plsd) / General (*.dat)'; ...
                                                '*.dat','General (*.dat)'; ...
                                                '*.*',  'All Files (*.*)'}, ...
                                                'Pick a file');
            
            
            if(~file_name)
                return;
            end

            full_name = [path_name, file_name];
            cd(path_name)

            [~,~,ext] = fileparts(full_name);

            output_folder_name = path_name;
            if(~output_folder_name)
                return;
            end
            ext = lower(ext);
            if strcmp(ext, '.edf')
                full_data_mat = edf2matlab2(full_name, output_folder_name, log);
            elseif strcmp(ext, '.txt')
                full_data_mat = etTxt2matlab(full_name, output_folder_name, log);
            elseif strcmp(ext, '.csv')
                full_data_mat = etCsv2matlab(full_name, output_folder_name, log);
            elseif strcmp(ext, '.asl')
                full_data_mat = asl2matlab(full_name, output_folder_name, log);
            elseif strcmp(ext, '.tbi')
                full_data_mat = tobii2matlab(full_name, output_folder_name, log);
            elseif strcmp(ext, '.plsd')
                full_data_mat = plsd2matlab(full_name, output_folder_name, log);

            elseif strcmp(ext, '.dat')
                full_data_mat = dat2matlab(full_name, output_folder_name, log);
            else
                full_data_mat = etTxt2matlab(full_name, output_folder_name, log);
            end
            if strfind(get(log, 'String'), 'Error')
                return;
            end
            data = read_data.parse_data(full_data_mat, log);
            print_log([strrep(file_name, '_', '\_') ' was loaded successfully! ' num2str(round(toc(start))) ' seconds. Sampling rate: ' num2str(data.rate) 'Hz'], log);    
        end

        function data = parse_data(full_data_mat, log)
            data = [];
            if(isempty(full_data_mat))
                return;
            end
            full_data_mat       = load_chap_data(full_data_mat, log);
            data.data_mat       = full_data_mat.data;
            vars                = data.data_mat.Properties.VariableNames();
            data.conds_data     = read_data.load_conditions(data.data_mat, vars, log);
            data.rate           = full_data_mat.rate;
            data.events2        = full_data_mat.events2;
            data.vars2          = full_data_mat.vars2;
            full_file_name = full_data_mat.file_name;
            [~, file_name, ext] = fileparts(full_file_name);
            data.file_name      = file_name;
            data.file_type      = ext(2:end); 
            
            full_data_mat.var_data_table.event_Trial_Onset = zeros(size(full_data_mat.var_data_table, 1), 1);
            data.var_data_table = full_data_mat.var_data_table;
        end
        
        function conds_data = load_conditions(data_mat, vars, log)
            buildin_var = {'pupil_x', 'pupil_y', 'pupil_size', 'trial_id', 'timestamps', 'msgs'};
            vars = vars(~ismember(vars, buildin_var));
            conds_data = [];
            index      = 0;
            print_log('Initiate condition''s levels loading', log);    
            
            for var=vars
                index    = index + 1;
                var_name = char(var);
                if(strfind(var_name, 'event_'))
                    vars = vars(~ismember(vars, var_name));
                    continue;
                end
                data_mat.(var_name)(cellfun(@isempty, data_mat.(var_name))) = {''};
                conds       = unique(data_mat.(var_name));
                conds       = conds(~strcmp(conds,''));
                rows        = max(size(conds,1), size(conds_data, 1));
                temp_arr    = cell(rows, size(conds_data, 2)+1);
                temp_arr(:) = {''};
                temp_arr(1:size(conds_data,1),1:size(conds_data,2)) = conds_data;
                temp_arr(1:size(conds,1),end) = conds;
                conds_data = temp_arr;
                percentage = round(100*(index/size(vars,2)));
                print_log(['Levels were loaded: ' num2str(percentage) '%'], log);    
            end
            print_log(['Loading of all variable levels has been completed: ' num2str(toc) ' seconds'], log);    
        end
    end
end