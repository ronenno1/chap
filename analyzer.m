classdef analyzer    
    methods(Static)
        % this function find all the posible combinations with the given conditions.
        % each number (0+) represents index of value in the j variable (in the NXM matrix).
        % for instance:
        % for given 3 variables X (ax, bx), Y (ay, by) Z (az, bz)
        % 0 1 1 - ay, az
        % 1 1 2 - ax, ay, bz
        % zeros mean unused variable. 
        function comps_ids = make_comps(conds)
            cond_names = fieldnames(conds);
            total_comps = 1;
            % create condition's names and calculate the number of them 
            for cond_name=1:size(cond_names, 1)
                conds.(char(cond_names(cond_name))) = vertcat('_',conds.(char(cond_names(cond_name))));
                total_comps = total_comps*size(conds.(char(cond_names(cond_name))),1); 
            end
            comps_ids      = zeros(total_comps,size(cond_names, 1));
            times_in_block = total_comps;

            for cond_name = 1:size(cond_names, 1)
                num_of_levels   = size(conds.(char(cond_names(cond_name))), 1);
                times_in_block  = times_in_block/num_of_levels;
                reps            = total_comps/(times_in_block*num_of_levels);
                block_size      = times_in_block*num_of_levels;
                for l=1:num_of_levels
                    comps_ids(1+(l-1)*times_in_block:l*times_in_block, cond_name) = l-1;
                end
                replicate_val = repmat(comps_ids(1:block_size,cond_name)',1,reps)';
                comps_ids(:, cond_name) = replicate_val;
            end
        end

        function data = analyze_all(data, data_mat, cond_names, comp_nums, conds, events_names, comp_names, log)
            counter = 0;
            for comp = 1:size(comp_nums, 1) % comp - index for posible comparation

                comp_name = [];
                used = true;
                for condition = 1:size(comp_nums, 2)
                    level     = comp_nums(comp, condition);
                    if(~level)
                        continue;
                    end
                    level_names = conds.(char(cond_names(condition)));
                    level_name  = level_names(level);
                    if(isempty(comp_name))
                        comp_name = level_name;
                    else
                        comp_name = strcat(comp_name, '_x_', level_name); 
                    end
                end
                comp_name = ['c_', char(strrep(strrep(char(comp_name), '.', '_'), '-', '_'))];
                if size(comp_names, 1)>0 && sum(strcmp(comp_names, comp_name))==0
                    continue;
                end

                comp_name = [];
                id = [];
                
                for condition = 1:size(comp_nums, 2)
                    cond_name = char(cond_names(condition));
                    level     = comp_nums(comp, condition);
                    if(~level)
                        continue;
                    end
                    level_names = conds.(char(cond_names(condition)));
                    level_name  = level_names(level);
                    if(isempty(comp_name))
                        comp_name = level_name;
                    else
                        comp_name = strcat(comp_name, '_x_', level_name); 
                    end
                    tmp_comp_name = ['c_', char(strrep(strrep(char(comp_name), '.', '_'), '-', '_'))];
                    if size(comp_names, 1)>0 && sum(not(cellfun('isempty',strfind(comp_names, tmp_comp_name))))==0
                        continue;
                    end
                    new_ids = find(strcmp(data_mat.(cond_name), level_name));
                    if(isempty(new_ids))
                        used = false;
                        continue;
                    end
                    if ~used
                        continue;
                    end

                    if(isempty(id)) && used
                        id = new_ids;
                    else
                        id = intersect(id, new_ids);
                    end
                    if(isempty(id))
                        used = false;
                        continue;
                    end

                end
                if ~used || isempty(id) 
                    continue;
                end
                
                comp_name = ['c_', char(strrep(strrep(char(comp_name), '.', '_'), '-', '_'))];
                counter   = counter + 1;

                total_comps = size(comp_nums, 1);
                if size(comp_names, 1)>0
                    total_comps = size(comp_names, 1);
                end
                [cond_mat_data, cond_events_data, cond_ids, avg_cond_mat, avg_events, outliers, valid_trials] = analyzer.parse_condition(data_mat, comp_name, events_names, id, counter, total_comps, log);

                data.cond_blinks.(comp_name)      = data.blinks_data(cond_ids);
                data.cond_mat_data.(comp_name)    = cond_mat_data;
                data.cond_ids.(comp_name)         = cond_ids;

                data.cond_events_data.(comp_name) = cond_events_data;
                data.cond_mat_cleaned.(comp_name) = avg_cond_mat;
                data.cond_events.(comp_name)      = avg_events;
                data.outliers.(comp_name)         = outliers;
                data.valid_trials.(comp_name)     = valid_trials;
            end    
            if (~exist('avg_cond_mat', 'var') || strcmp(log, ''))
                return;
            end
            print_log('Parsing of all levels has been completed' , log);    

        end
        
        function [cond_mat_data, cond_events_data, trials, avg_cond_mat, avg_events, outliers, valid_trials] = parse_condition(data_mat, comp_name, events_names, id, counter, total_comps, log)
            cond_mat.(comp_name)                = data_mat(id, :); 
            [trials, trial_ids]                 = unique(cond_mat.(comp_name).trial_id);

            trials_data.(comp_name).trials      = trials;
            trials_data.(comp_name).trial_ids   = trial_ids;
            [cond_mat_data, cond_events_data, avg_cond_mat, avg_events, outliers, valid_trials]  = analyzer.get_condition_avgs(cond_mat.(char(comp_name)), trials_data.(char(comp_name)), events_names);
            print_log(['Levels were parsed: ' num2str(round(100*(counter/total_comps))) '%'], log);    
        end

        function [cond_mat_data, cond_events_data, avg_cond_mat, avg_events, outliers, valid_trials] = get_condition_avgs(cond_mat, trials_data, event_names)
            avg_events  = [];
            trial_ids   = trials_data.trial_ids;
            trials      = trials_data.trials;
            %% get the max number of samples (the longest trial)
            max_samples = max(diff(trial_ids));
            for evnt = 1:size(event_names, 2)
                avg_events.(char(event_names(:,evnt))) = 0;
            end
            cond_mat_data = zeros(max_samples, size(trials, 1));
            %% get data
            Trial_Offset   = 0;
            outliers    = 0;
            valid_trials = 0;
            num_of_trials = size(trials, 1);
            
            cond_events_data = zeros(size(event_names, 2)+1, num_of_trials);
            for trial=1:num_of_trials
                %% get relevant samples
                if trial < size(trials, 1)
                    trial_data = cond_mat(trial_ids(trial):trial_ids(trial+1)-1, :);
                else
                    trial_data = cond_mat(trial_ids(trial):end, :);
                end

                %% get values for events
                Trial_Offset = Trial_Offset + mean(trial_data.event_Trial_Offset);
                for evnt = 1:size(event_names, 2)
                    if (evnt==1)
                    end
                    cond_events_data(evnt, trial) = cond_mat.(char(event_names(:,evnt)))(trial_ids(trial));
                    avg_events.(char(event_names(:,evnt))) = avg_events.(char(event_names(:,evnt))) + cond_mat.(char(event_names(:,evnt)))(trial_ids(trial));
                end
                cond_events_data(end, trial) = cond_mat.event_Trial_Offset(trial_ids(trial));

                %% find outliers
                if(max(trial_data.pupil) > 0)
                    valid_trials = valid_trials+1;
                else
                    outliers = outliers + 1;
                end
                cond_mat_data(1:size(trial_data.pupil, 1), trial) = trial_data.pupil;
            end
            cond_mat_data(cond_mat_data==0)=nan; 
            for evnt = 1:size(event_names, 2)
                avg_events.(char(event_names(:,evnt))) = avg_events.(char(event_names(:,evnt)))/num_of_trials;
            end
            
            avg_events.Trial_Offset   = Trial_Offset/num_of_trials;
            avg_cond_mat = mean(cond_mat_data, 2, 'omitnan');
        end
    end
end

