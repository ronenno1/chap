classdef ploter2   
    methods(Static)
        function [ploted_data, events] = draw_graph(cond_mat_data, cond_events_data, blinks_data, comp_names, rate, fig, data_mean, data_std, configuration, file_name)
            if (~exist('file_name', 'var'))
                file_name = '';
            end
            if (~isfield(configuration, 'relative_val'))
                configuration.from_val               = '';
                configuration.to_val                 = '';
                configuration.BinsNumber_val         = 0;
                configuration.relative_val           = 'no';
                configuration.Method_val             = 'arbitrary units';
                configuration.PreEventNumber_val     = 0;
                configuration.baseline_val           = 0;
                configuration.scattering_val         = 'no';

            end

%             if (~strcmp(configuration.relative_val, 'percentage'))
%                 configuration.Method_val = 'arbitrary units';
%             end
            [ploted_data, events] = ploter2.do_draw_graph(cond_mat_data, cond_events_data, blinks_data, comp_names, configuration.event_names, rate, fig, data_mean, data_std, configuration.BinsNumber_val, configuration.PreEventNumber_val, configuration.Method_val, configuration.from_val, configuration.to_val, configuration.relative_val, configuration.baseline_val, configuration.scattering_val, file_name);
        end

        % for each condition we will take the cleaned data and then:
        % 1.  cut it from the begining (the time before the first event) until the first event.
        % 2.  convert it to the relevant method (mm or Z-score)
        % 3.  do the same (1+2) from the first event until the last event.
        % 4.  (optional) convert the data before (and after) the first event to bins.
        % 5.  define the value in the first event point to be the avarage between the last point before the event and the first point after the event.
        % 6.  create the x-axis according to the relevant range.
        % 7.  (optional) calculate relative values for each point to the first event.
        % 8.  concat the values before and after the first event.
        % 9.  draw the curve with color number n (according to the number of the condition)
        % 10. save the curve and the x-axis values to ploted_data.
        % 11. draw legend.
        function [ploted_data, events_id_to_show, valid_comps] = do_draw_conditions(cond_mat_data, cond_events_data, blinks_data, comp_names, event_names, cmap, ms, fig, data_mean, data_std, bins_ms, PreEventNumber, method, from, to, relative, baseline, scattering)
            %debug
%             from = 'Frame_to_be_displayed_1';
%             to = 'event_Frame_to_be_displayed_236';
%             relative = 'percentage';
%             method = 'Z-score';
            Bins = round(max(bins_ms, 1)/ms);
            if(Bins == Inf || Bins == 0 )
                Bins=1;
            end
            
            pre_event_samples = round(floor(PreEventNumber/ms)/Bins);
            
            if(baseline == 0)
                baseline = 1;
            end

            baseline_samples   = max(1, round(baseline/(ms*Bins)));
            
            from_event_id = [];
            to_event_id = [];
            if(~strcmp(char(from),'')&&~strcmp(char(from),'Trial_Onset'))
                from_event_id = find(ismember(event_names, strcat('event_', char(from))));
            end
            if(~strcmp(char(to),'')&&~strcmp(char(to),'Trial_Offset'))
                to_event_id   = find(ismember(event_names, strcat('event_', char(to))));
            end
            
            if(isempty(from_event_id))
                first_event_to_show = 0;
            else
                first_event_to_show = from_event_id;
            end
            
            last_event_to_show  = to_event_id;
            if(isempty(last_event_to_show))
                last_event_to_show = size(event_names, 2)-1;
            end
            if strcmp(event_names, 'Trial_Offset')
                last_event_to_show = 1;
            end
            
            events_id_to_show = first_event_to_show:last_event_to_show;
            
            cmap2 = [0, 255, 0;
                    255, 0, 0; 
                    0, 0, 255;
                    255, 0, 255;
                    0, 255, 255;
                    255, 255, 0; 
                    0, 128, 0;
                    128, 0, 0; 
                    0, 0, 128;
                    128, 0, 128;
                    0, 128, 128;
                    128, 128, 0; 
                    128, 128, 128; 
                    0, 0, 0]/255;
                
            cmap3 = cmap2+0.8;
            cmap3(cmap3>1) = 1;
            hold on;
            valid_comp_id = 0;
            valid_comps = {''};
            for comp = 1:size(comp_names, 1)

                cond_name = char(comp_names(comp));
                if~isfield(cond_mat_data, cond_name)
                    continue;
                end
                if ~isnan(mean(mean(cond_mat_data.(cond_name), 'omitnan'), 'omitnan')) 
                    valid_comp_id = valid_comp_id +1;
                    valid_comps{valid_comp_id, :} = cond_name;
                end

                % cut, calculate the baseline and align to the first event
                cuted_data = zeros(round(size(cond_mat_data.(cond_name), 1)/Bins), size(cond_mat_data.(cond_name), 2));

                for trial=1:size(cond_mat_data.(cond_name), 2)
                    trial_full_data = cond_mat_data.(cond_name)(:, trial);
                    trial_full_data = converter.convert2method(trial_full_data, method, data_mean, data_std);
                    if(Bins>1)
                        bins_data  = arrayfun(@(x) mean(trial_full_data(x:min(x+Bins-1, end), :), 1), 1:Bins:size(trial_full_data, 1), 'UniformOutput', false);
                        trial_full_data = cat(1, bins_data{:});
                    end
                    start_from = 1;
                    if(~isempty(from_event_id))
                        start_from = round(cond_events_data.(cond_name)(from_event_id, trial)/Bins); %%TODO add ms
                    end
                    end_with = size(trial_full_data, 1);
                    

                    if(~isempty(to_event_id))
                        end_with = round(cond_events_data.(cond_name)(to_event_id, trial)/Bins);
                    end
                    
                    
                    if(start_from<=0)
                        continue;
                    end
                    baseline_data = mean(trial_full_data(max(1, start_from-baseline_samples+1):start_from), 'omitnan');
                    if(start_from && start_from-pre_event_samples<=0)
                        continue;
                    end
                    
                    trial_data = trial_full_data(max(1, start_from-pre_event_samples):end_with);
                    if(~strcmp(relative, 'no'))
                        if(strcmp(relative, 'percentage'))
                            trial_data = 100*trial_data/baseline_data-100;
                        else
                            trial_data = trial_data-baseline_data;
                        end
                    end
                    cuted_data(1:end_with-start_from+1+pre_event_samples, trial) = trial_data;                    
                end
                
                events      = cond_events_data.(cond_name);
                
                evant_avg = mean(events,2);
                end_with_avg =  round(floor(evant_avg(end))/Bins);

                if(~isempty(to_event_id))
                    end_with_avg =  round(floor(evant_avg(to_event_id))/Bins);
                end

                start_with_avg =  1;

                if(~isempty(from_event_id))
                    start_with_avg =  round(floor(evant_avg(from_event_id))/Bins);
                end
                
                % avaraging all the trials
                cuted_data(pre_event_samples+1) = cuted_data(pre_event_samples+1) + 0.001;
                cuted_data(cuted_data==0)=nan;
                cuted_data(pre_event_samples+1) = cuted_data(pre_event_samples+1) - 0.001;
                avg_cond_mat = mean(cuted_data, 2, 'omitnan');
                std_cond_mat = std(cuted_data', 'omitnan')';
                num_of_trials = size(cuted_data, 2);

                % remove this condition in case of one trial
%                 if(num_of_trials <2)
%                     continue;
%                 end
                
                if(end_with_avg-start_with_avg+pre_event_samples>length(cuted_data))
                    continue;
                end
                if(start_with_avg<pre_event_samples)
                    continue;
                end

                data4analyze = cuted_data(1:end_with_avg-start_with_avg+pre_event_samples,:);
                avg_cond_mat = avg_cond_mat(1:end_with_avg-start_with_avg+pre_event_samples,:);  
                
                if num_of_trials<2
                    std_cond_mat = zeros(size(data4analyze));                    
                else
                    std_cond_mat = std_cond_mat(1:end_with_avg-start_with_avg+pre_event_samples,:);                    
                end
                

                
                avg_cond_mat(find(isnan(avg_cond_mat))) = [];
                std_cond_mat(find(isnan(std_cond_mat))) = [];
                % convert the data to the relevant units

                % create x_axis
                Trial_Offset_ms = round(size(avg_cond_mat, 1)*ms*Bins);
                x_axis = linspace(0 ,Trial_Offset_ms, size(avg_cond_mat,1))-PreEventNumber;
               
                
                x = x_axis';
                
                
                d  = avg_cond_mat;
                sd = std_cond_mat;
                n  = [];
                ci = sd;
                se = sd;
                for i=1:size(data4analyze, 1)
                    n(i) = sum(~isnan(data4analyze(i, :)));
                end
                if n ~=0
                    se = sd./(n.^0.5)';
                    for i=1:size(data4analyze, 1)
                        ci(i, :) = tinv(0.975, n(i)-1)*se(i);
                    end
                end
%                 ci = tinv(0.975, n-1)*se;
                ci_data = cuted_data(1:end_with_avg-start_with_avg+pre_event_samples,:);                    

                q = 0.95;
                HDImin = zeros(length(x), 1);
                HDImax = zeros(length(x), 1);

                for s = 1:size(ci_data,1)
                    try
                        samples = ci_data(s, :);
                        samples(find(isnan(samples))) = [];
                        sortedVec = sort(samples);
                        ciIdx = ceil(q * length(sortedVec));
                        nCIs = length(sortedVec) - ciIdx;  % number of vector elements that make HDI

                        % Determine middle of HDI to get upper and lower bound
                        ciWidth = zeros(nCIs,1);
                        for ind = 1:nCIs
                            ciWidth(ind) = sortedVec(ind + ciIdx) - sortedVec(ind);
                        end

                        [~,idxMin] = min(ciWidth);
                        HDImin(s, :) = sortedVec(idxMin);
                        HDImax(s, :) = sortedVec(idxMin + ciIdx);
                    catch
                        HDImin(s, :) = nan;
                        HDImax(s, :) = nan;
                    end
                end
                
%                 h = fill([x;flipud(x)],[HDImin; flipud(HDImax)], cmap3(comp,:), 'EdgeColor', cmap2(comp,:));

                ploted_data.(cond_name).pupil      = avg_cond_mat;
                ploted_data.(cond_name).x_axis     = x_axis;
                ploted_data.(cond_name).blinks     = blinks_data.(cond_name);
                ploted_data.(cond_name).cuted_data = cuted_data;

                hold on
                if(size(cond_events_data.(char(cond_name)),1)==0)
                    continue;
                end
                if size(x_axis, 2)>1
                    xlim(fig, [x_axis(1) x_axis(end)])
                end

                if strcmp(scattering, 'no')
                   continue; 
                end
                
                if strcmp(scattering, 'SD')
                    scattering_data = sd;
                elseif strcmp(scattering, 'SE')
                    scattering_data = se;
                elseif strcmp(scattering, 'CI')
                    scattering_data = ci;
                end
%                 h = fill([x;flipud(x)],[d-scattering_data;flipud(d+scattering_data)], cmap3(comp,:), 'LineStyle', '--', 'EdgeColor', cmap2(comp,:));
%                 set(h,'facealpha',.2)
                h = fill([x;flipud(x)],[d-scattering_data;flipud(d+scattering_data)], cmap3(comp,:), 'LineStyle', '-', 'EdgeColor', cmap2(comp,:), 'Parent', fig);
                set(h,'facealpha',.5)
                
            end
            if (~exist('ploted_data', 'var'))
                ploted_data.(cond_name).pupil      = [];
                ploted_data.(cond_name).x_axis     = [];
                ploted_data.(cond_name).blinks     = [];
                ploted_data.(cond_name).cuted_data = [];
                return;
            end    
            %origin
            do_diff = 0;
            if do_diff
                for comp = 1:size(comp_names, 1)/2
                    cond_name1 = char(comp_names(comp));
                    cond_name2 = char(comp_names(comp+2));
                    min_dur = min(length(ploted_data.(cond_name1).pupil), length(ploted_data.(cond_name2).pupil));
                    diff_pupil = ploted_data.(cond_name2).pupil(1:min_dur)-ploted_data.(cond_name1).pupil(1:min_dur);
                    ploted_data.(cond_name1).pupil = diff_pupil;
                    ploted_data.(cond_name1).x_axis = ploted_data.(cond_name1).x_axis(1:min_dur);
                    ploted_data.(cond_name2).x_axis = ploted_data.(cond_name1).x_axis(1:min_dur);
                    ploted_data.(cond_name2).pupil = diff_pupil*0;
                end
            end

            for comp = 1:size(comp_names, 1)
                cond_name = char(comp_names(comp));
                if~isfield(ploted_data, cond_name)
                    continue;
                end
                plot(ploted_data.(cond_name).x_axis, ploted_data.(cond_name).pupil, 'LineWidth',3, 'Color', cmap2(comp,:), 'Parent', fig);
            end            
            set(fig, 'LineWidth', 2);
            try
                xtickformat(fig, '%,.4g');
            catch
            end
        end
        
        function [ploted_data, events] = do_draw_events(ploted_data, cond_events_data, comp_names, event_names, cmap, ms, events_id_to_show, fig)
            cmap2 = [0, 255, 0;
                    255, 0, 0; 
                    0, 0, 255;
                    255, 0, 255;
                    0, 255, 255;
                    255, 255, 0; 
                    0, 128, 0;
                    128, 0, 0; 
                    0, 0, 128;
                    128, 0, 128;
                    0, 128, 128;
                    128, 128, 0; 
                    128, 128, 128; 
                    0, 0, 0]/255;

            
            relative = true;
            first_event = events_id_to_show(1);
            if(first_event<1)
                events_id_to_show = events_id_to_show(2:end);
                first_event = events_id_to_show(1);
                relative = false;
            end
            if ~iscell(event_names)
                event_names = cellstr(event_names);
            end
            for comp = 1:size(comp_names, 1)
                cond_name = char(comp_names(comp));
                if~isfield(cond_events_data, cond_name)
                    continue;
                end
                if(size(cond_events_data.(char(cond_name)), 1)==0)
                    continue;
                end
                events    = cond_events_data.(cond_name);
                all_events    = events*ms;
                outlier_trials = find(isnan(min(ploted_data.(char(cond_name)).cuted_data)));
                events(:, outlier_trials) = [];
                ploted_data.(char(cond_name)).outlier_trials = outlier_trials;

                evant_avg = mean(events, 2);
                cond_events = events;
                cond_events_ms = ms*cond_events;
                for e = events_id_to_show  
                    ploted_data.(cond_name).all_events(e, :) = events(e, :);
                    ploted_data.(cond_name).all_events_ms(e, :) = events(e, :)*ms;


                    event_time = evant_avg(e);
                    if(relative)
                        event_time = event_time-evant_avg(first_event);
                    end
                    y = get(fig, 'ylim');
                    event_time_ms = ms*event_time;

                    ploted_data.(cond_name).events.(char(event_names(e))) = event_time_ms;
                    color2plot = [0, 0, 0];
                    if event_time_ms~=0
                        color2plot = cmap2(comp,:);
                    end
                    plot([event_time_ms event_time_ms], y, 'Color', color2plot, 'LineStyle', '--', 'Parent', fig);
                end
                event_time = evant_avg(end);
                if(relative)
                    event_time = event_time-evant_avg(first_event);
                end
                event_time_ms = ms*event_time;
                ploted_data.(cond_name).events.Trial_Offset = event_time_ms;
            end
        end
        
        function do_draw_lables(relative, Method, file_name, fig)
            y_unit_label        = 'arbitrary units';
            relative_label      = '';
            if strcmp(relative, 'difference') || strcmp(relative, 'percentage')
                relative_label = 'Relative ';
            end
            if strcmp(Method,'mm')
                y_unit_label = 'mm';
            end
            if strcmp(Method,'Z-score')
                y_unit_label = 'Z-score';
            end
            if strcmp(relative, 'percentage')
                y_unit_label = '% ';
            end
            xlabel(fig, 'Time [ms]', 'FontWeight','bold');
            ylabel(fig, [relative_label 'Pupil Size [' y_unit_label ']'], 'FontWeight','bold');
            set(gca,'FontWeight','bold');            
            set(gca,'box','on');

            title(fig, strrep(file_name, '_', '\_'));        
        end
        
        function [ploted_data, events] = do_draw_graph(cond_mat_data, cond_events_data, blinks_data, comp_names, event_names, rate, fig, data_mean, data_std, Bins, PreEventNumber, Method, from, to, relative, baseline, scattering, file_name)
            cla(fig);
            cmap = hsv(size(comp_names, 1));
            ms   = 1000/rate;
            [ploted_data, events_id_to_show, valid_comps] = ploter2.do_draw_conditions(cond_mat_data, cond_events_data, blinks_data, comp_names, event_names, cmap, ms, fig, data_mean, data_std, Bins, PreEventNumber, Method, from, to, relative, baseline, scattering);
            
            if length(fieldnames(ploted_data))<length(comp_names)
                events = [];
                ploted_data.wrong_range = true;
                return
            end

            [ploted_data, events] = ploter2.do_draw_events(ploted_data, cond_events_data, comp_names, event_names, cmap, ms, events_id_to_show, fig);
            ploter2.do_draw_lables(relative, Method, file_name, fig);
            
            comp_names_fixed = cellfun(@(x) x(3:end), valid_comps, 'UniformOutput', false);
            comp_names_fixed = strrep(strrep(comp_names_fixed, '_x_', ' & '),'_',' ');
            legend(fig, char(comp_names_fixed), 'Location', 'Best');
            if ~strcmp(scattering, 'no')
                chs = flipud(get(fig,'children'));
                legend(chs(size(comp_names_fixed, 1)+1:size(comp_names_fixed, 1)*2), char(comp_names_fixed), 'Location', 'Best');
            end

        end
    end
end
