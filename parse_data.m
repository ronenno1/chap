classdef parse_data
    methods(Static)
        function var_data_table = parse_external_vars(vars_file)
            if exist('detectImportOptions', 'file')
                var_data_table = readtable(vars_file, detectImportOptions(vars_file));
            else
                var_data_table = readtable(vars_file, 'Delimiter', ',');
            end

            
            if ~ismember('trial_id', var_data_table.Properties.VariableNames)
                error('Corrupt variable file (trial\_id column does not exist)');
            end
            var_data_table.trial_id = [];
        end
        
        function total_var_data = parse_vars(event_msgs, event_timestamps, timestamps, var_ids, Trial_Onset_num, external)
            if (~exist('external', 'var'))
                external = false;
            end
            var_full_data = event_msgs(var_ids);
            for var = 1:size(var_full_data, 1)
                timestamp = event_timestamps(var_ids(var));
                sample_id = find(timestamps<=timestamp, 1, 'last');
                trial_id  = find(Trial_Onset_num<sample_id, 1, 'last');
                var_data_arr  = strsplit(cell2mat(var_full_data(var)));
                if external
                    total_var_data.(char(var_data_arr(1)))(trial_id, :) = var_data_arr(2);
                else
                    if size(var_data_arr, 2)==3
                        var_data_arr{4} = '-';
                    end
                    total_var_data.(regexprep(char(var_data_arr(3)), '\W', '_'))(trial_id, :) = strrep(var_data_arr(4), '''', '') ;
                end
            end
            var_names = fieldnames(total_var_data);
            max_values = size(Trial_Onset_num, 1);

            for var = 1:size(var_names, 1)
                if size(total_var_data.(char(var_names(var))), 1) < max_values
                    total_var_data.(char(var_names(var))){max_values, :} = [];
                end
            end
        end
        
        function event_data_table = parse_events(event_msgs, event_timestamps, timestamps, Trial_Onset_num, Trial_Offset_num, external)
            if (~exist('external', 'var'))
                external = false;
            end

            for event=1:size(event_msgs, 1)
                timestamp  = event_timestamps(event);
                sample_id  = find(timestamps>=timestamp, 1, 'first');
                trial_id   = find(Trial_Onset_num<=sample_id, 1, 'last'); 

                event_str  = char(event_msgs(event));
                if(external)
                    event_data = strsplit(event_str);
                    event_name  = ['event_' char(event_data(2:end))];
                else
                    event_data = strsplit(event_str);
                    event_name = ['event_' char(event_data(3))];
                end
                if(isempty(trial_id))
                    event_parsed_data.trial_id   = -1;
                    event_parsed_data.event_val  = -1;
                    continue;
                end
                if external && length(event_data) == 2 && str2double(event_data{1})>0
                    event_diff_val = (timestamps(Trial_Offset_num(trial_id))-timestamps(Trial_Onset_num(trial_id))-str2double(event_data{1}));
                else
                    event_diff_val = (timestamp-timestamps(Trial_Onset_num(trial_id)));
                end
                total_event_data.(regexprep(char(event_name), '\W', '_'))(trial_id, :) = event_diff_val;    
            end
            event_names = fieldnames(total_event_data);
            max_values = size(Trial_Onset_num, 1);
            for event = 1:size(event_names, 1)
                if size(total_event_data.(char(event_names(event))), 1) < max_values
                    total_event_data.(char(event_names(event)))(max_values, :) = 0;
                end
            end

            event_data_table = struct2table(total_event_data);
        end
   
        function [timestamps, data] = add_missing_samples(timestamps, data, log)
            diffs = diff(timestamps);
            misses = find(diffs>1/data.rate);
            % if there are too many missing samples, it is highly recommended to
            % add 0.00001 to the sampling rate:
            misses = find(diffs>1/data.rate+0.00001);
            gap_counter = 0;
            percentages = zeros(1, 10);

            for miss = 1: length(misses)
                percentage = round(100*(miss/length(misses)));
                
                if percentage>0 && ~mod(percentage, 10) && percentages(percentage/10)==0 
                    print_log(['Sampales were added: ' num2str(percentage) '%'], log);    
                    percentages(percentage/10) = 1;
                end

                if misses(miss)+2>length(timestamps)
                    continue;
                end

                gaps2 = round((timestamps(misses(miss)+2)-timestamps(misses(miss)))*data.rate)-2;
                gaps1 = round((timestamps(misses(miss)+1)-timestamps(misses(miss)))*data.rate)-1;

                gaps  = min(gaps1, gaps2);
                if miss<length(misses) && misses(miss+1)==misses(miss)+1
                    gaps = round((timestamps(misses(miss)+1)-timestamps(misses(miss)))*data.rate)-1;
                end
                if gaps<0
                    continue;
                end

                gap_id = gaps;
                gap_counter = gap_counter+gap_id;
                while gap_id 

                    [gap_id, misses(miss)];
                    trial_time   = timestamps(misses(miss))+(1/data.rate);

                    timestamps = [timestamps(1:misses(miss)), trial_time, timestamps(misses(miss)+1:end)]; 
                    new_time = datestr(((trial_time/86400)+ datenum('01-Jan-1970')), 'yyyy-mm-dd HH:MM:SS.FFF');

                    data.timestamps = [data.timestamps(1:misses(miss)); new_time; data.timestamps(misses(miss)+1:end)]; 

                    data.pupil_size = [data.pupil_size(1:misses(miss)); 0; data.pupil_size(misses(miss)+1:end)]; 
                    data.pupil_x = [data.pupil_x(1:misses(miss)); 0; data.pupil_x(misses(miss)+1:end)]; 
                    data.pupil_y = [data.pupil_y(1:misses(miss)); 0; data.pupil_y(misses(miss)+1:end)]; 
                    gap_id = gap_id-1;
                    misses(miss:end)=misses(miss:end)+1; 
                end        
            end
            data.gap_counter = gap_counter/data.rate;

            %% recheck
%             diffs = diff(timestamps);
%             missesa = find(diffs>=1/data.rate);
% 
%             gap_counter = 0;
%             for miss = 1: length(missesa)
%                 if misses(miss)+2>length(timestamps)
%                     continue;
%                 end
% 
%                 gaps = round((timestamps(misses(miss)+2)-timestamps(misses(miss)))*data.rate)-2;
% 
%                 if miss<length(misses) && misses(miss+1)==misses(miss)+1
%                     gaps = round((timestamps(misses(miss)+1)-timestamps(misses(miss)))*data.rate)-1;
%                 end
% 
%                 gap_counter = gap_counter+gaps;
%             end
        end
    end
end