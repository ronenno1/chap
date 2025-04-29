classdef blinks_analysis
    methods(Static)        
        
        function run_blinks_analysis()
            [data, folder, files, condition2analyze, condition_names] = blinks_analysis.load_files();
            if isempty(data)
                return;
            end
            close(findobj('type', 'figure', 'name', 'CHAP - Eye-blinks Analysis'))
        
            figure = gui_lib.create_figure('CHAP - Eye-blinks Analysis', 'images/chap_3.jpg', [0 0 900 750]); %create gui figure
            gui_lib.uicontrol_title('CONFIGURATION', 130, 70);
            gui_lib.uicontrol_title('RESULTS', 578, 70);
            gui_lib.uicontrol_title('RESULTS', 418, 493);
            log   = text(35, 730, '', 'HorizontalAlignment', 'left', 'Color', [1 1 1], 'FontWeight', 'bold');
            
            
            hp	  = uipanel('Parent',figure, 'Units','Pixels', 'BackgroundColor',[1 1 1], 'Position',[362 325 500 340]);    
            fig   = subplot(1, 1, 1, 'Parent', hp);
            analyze_bot = gui_lib.uicontrol_button(figure, [32 281 300 30], 'Analyze');
    

            exit_bot = gui_lib.uicontrol_button(figure, [32 40 125 30], 'Exit');
            set(exit_bot, 'callback', {@(~,~)blinks_analysis.exit});
            
            
            tooltip = '<html><i>Save data';
            gui_lib.uicontrol_button(figure, [162 40 155 30], 'Save data', {@(~,~)blinks_analysis.save_stat(figure, folder)}, tooltip); 

            
            data = blinks_analysis.show_statistical_vars(data, figure);
            data.folder = folder;
            
            axis = data.x_axis;
            set(log, 'String', '');
            if isempty(files)
                print_log('Error: Files cannot be found!', log);    
                return;
            end
            
            print_log('Parsing data...', log);    

            blinks = blinks_analysis.parse_all_data([folder filesep, 'mat'], files, condition2analyze, log);

            print_log('Saving data table...', log);    

            blinks_analysis.save_data_table(blinks(:, 2), condition_names, axis, folder, files);

            set(analyze_bot, 'callback', {@blinks_analysis.calc_stat_params, data, blinks(:, 2), figure});
            
            contrasts = dir([folder, 'stat', filesep '*vs*.csv']);
            contrasts = {contrasts.name}';

            contrasts = cellfun(@(x) x(1:end-4), contrasts, 'UniformOutput', false); 
            
            hold on;
            
            print_log('Printing data...', log);    
        
            blinks_analysis.print_figure(blinks(:, 2), condition_names, axis, fig);
            print_log('Calculating contrasts...', log);    

            blinks_analysis.print_contrasts(blinks, contrasts, data.configuration.comp_names, axis, folder, fig, log);
            
            blinks_analysis.calc_stat_params(0, 0, data, blinks(:, 2), figure);
            ylims = get(fig, 'ylim');
            file_name = [folder, 'EBR', filesep, 'EBR'];

            output.save_figure2(fig, [file_name, '.fig'], ylims);
            output.save_figure2(fig, [file_name, '.png'], ylims);
            print_log('Done!', log);    

        end
        
        function save_data_table(blinks, condition_names, axis, folder, files)
            path2save = [folder, filesep, 'EBR'];
            if ~exist([folder, filesep, 'EBR'], 'file')
                mkdir(path2save);
            end
            bin = zeros(length(axis), length(condition_names));
            for cond_id = 1:length(condition_names)
                bin = blinks{cond_id};
                files     = cellfun(@(x) x(1:end-4), files, 'UniformOutput', false);
                
                bin = bin(:, 1:length(axis));
                bin = [axis; bin];
                blinks_data  = array2table(bin);
                blinks_data = [array2table(['time'; files], 'VariableNames', {'Id'}), blinks_data];
                writetable(blinks_data, [path2save, filesep, condition_names{cond_id}, '.csv']);
            end
        end
        
        function [data, base_path, files, condition2analyze, condition_names] = load_files()
            data              = [];
            files             = [];

            condition2analyze = [];
            condition_names   = [];
            [file_name, base_path]  = uigetfile('*.mat', 'Select mat file ');
            if ~base_path
               disp('Error: *.mat file have to be chosen');
               return;
            end

            full_mat_name = [base_path, file_name];

            data_file = load(full_mat_name);
            if ~isfield(data_file, 'total_data')
               disp('Error: Wrong mat file');
               return;
            end
            data = data_file.total_data;
            condition2analyze = data.configuration.comp_names;

            condition_names = cellfun(@(x) x(3:end), condition2analyze, 'UniformOutput', false);
            condition_names = strrep(strrep(condition_names, '_x_', ' & '), '_', ' ');
            
            folder = [base_path, filesep, 'mat'];
            csv_data = readtable([base_path, filesep, 'csv', filesep, 'trials_data.csv']);
            first_event = strcat('event_', data.configuration.from_val);
            data.mean_first_event = round(mean(csv_data.(char(first_event))));
            data.events_data = csv_data;

            files = dir([folder, filesep, '*.mat']);
            files = {files.name}';
           
        end

        function blinks = parse_all_data(folder, files, condition2analyze, log)
            tic
            blinks = cell(length(condition2analyze), 2);
            cond_names = cellfun(@(x) x(3:end), condition2analyze, 'UniformOutput', false);
            cond_names = strrep(strrep(cond_names, '_x_', ' & '),'_',' ');

            for condition = 1:length(condition2analyze)
               max_duration = 0;
                for file = 1:length(files)
                    data_file    = load([folder, filesep, files{file}]);
                    cond_data    = data_file.ploted_data.(condition2analyze{condition}).cuted_data';
                    max_duration = max(max_duration, size(cond_data, 2));
                end
            end
            file_id = 0;
            percentages = zeros(1, 10);
            for condition = 1:length(condition2analyze) 
                mean_cond = zeros(file, max_duration);

                for file = 1:length(files)
                    file_id = file_id+1;
                    
                    percentage = round(100*(file_id /(length(files)*length(condition2analyze))));
                    if percentage>0 && ~mod(percentage, 10) && percentages(percentage/10)==0 
                        print_log(['Parsing data files: ' num2str(percentage) '%'], log);    
                        percentages(percentage/10) = 1;
                    end

                    data_file = load([folder, filesep, files{file}]);

                    cond_data = data_file.ploted_data.(condition2analyze{condition}).blinks;
                    outliers = data_file.ploted_data.(condition2analyze{condition}).outlier_trials;
                    cond_data(outliers) = [];
                    cond_blinks = zeros(length(cond_data), max_duration);

                    x_axis = data_file.ploted_data.(condition2analyze{condition}).x_axis;
                    for trial = 1:size(cond_blinks, 1)
                        trial_blinks = cond_data{trial};
                        blink_id = 1;
                        
                        while blink_id+1<=length(trial_blinks)
                            last_sample = min(1+trial_blinks(blink_id+1), max_duration);
                            cond_blinks(trial, 1+trial_blinks(blink_id): last_sample) = 1;
                            blink_id = blink_id + 2;
                        end

                        first_sample_ms = round(data_file.ploted_data.(condition2analyze{condition}).all_events(:, trial));
                        first_sample_ms = first_sample_ms(1);
                        first_sample = round((first_sample_ms-(x_axis(1)))/(1000/data_file.ploted_data.rate));
                        cond_blinks(trial, :) = [cond_blinks(trial, first_sample:end), zeros(1, first_sample-1)];
                    end
                    mean_cond(file, :) = mean(cond_blinks);
                end
                blinks{condition, 1} = cond_names{condition};
                blinks{condition, 2} = mean_cond;
            end
            print_log(['Data has been parsed successfully! ', num2str(toc) ' seconds'], log);    
%             diff = 0;
%             if diff
%                for cond = 1:size(blinks, 1)/2
%                    blinks{cond, 2} = blinks{cond, 2}- blinks{cond+2, 2};
%                    blinks{cond+2, 2} = blinks{cond, 2}*0;
%                end
%            end
        end

        
        function print_contrasts(blinks, contrasts, comp_names, axis, folder, fig, log)
            tic
            ylims =   get(gca,'ylim');
            min_val = ylims(1);
            max_val = ylims(2);

            range = max_val-min_val;
            num_of_contrasts = length(contrasts);
            sample_id = 0;
            percentages = zeros(1, 10);
%             files = dir([folder, 'EBR', filesep '*.csv']);
%             files = {files.name}';
%             data_files =cell(length(files), 2);
%             for file = 1:length(files)
%                 data_files{file, 1} = files{file};
%                 data_files{file, 2} = readtable([folder, 'EBR', filesep, files{file}]);
%             end
            for contranst_id   = 1:length(contrasts)
                contrast_arr   = split(contrasts{contranst_id}, ' vs ');
%                 first_cond     = readtable([folder, 'EBR', filesep, contrast_arr{1}, '.csv']);
%                 first_cond.Id  = [];
%                 first_cond     = table2array(first_cond);
%                 first_cond     = first_cond(2:end, :);
                first_cond_id = find(cellfun(@(x) strcmp(x, contrast_arr{1}), blinks(:, 1), 'UniformOutput', 1));
                first_cond    = blinks{first_cond_id, 2}; 
                first_cond    = first_cond(:, 1:length(axis));
%                 second_cond    = readtable([folder, 'EBR', filesep, contrast_arr{2}, '.csv']);
%                 second_cond.Id = [];
%                 second_cond    = table2array(second_cond);
%                 second_cond    = second_cond(2:end, :);
                second_cond_id = find(cellfun(@(x) strcmp(x, contrast_arr{2}), blinks(:, 1), 'UniformOutput', 1));
                second_cond    = blinks{second_cond_id, 2}; 
                second_cond    = second_cond(:, 1:length(axis));
                stat_data.time     = zeros(length(axis), 1);
                stat_data.tValues = zeros(length(axis), 1);
                stat_data.pValues = zeros(length(axis), 1);
                stat_data.BFs     = zeros(length(axis), 1);
                stat_data.pes     = zeros(length(axis), 1);
                stat_data.sd      = zeros(length(axis), 1);

                for sample = 1:length(axis)
                    sample_id = sample_id+1;
                    percentage = round(100*(sample_id /(length(axis)*length(contrasts))));
                    if percentage>0 && ~mod(percentage, 10) && percentages(percentage/10)==0 
                        print_log(['Analyzing data: ' num2str(percentage) '%'], log);    
                        percentages(percentage/10) = 1;
                    end

                    first_sample  = first_cond(:, sample);
                    second_sample = second_cond(:, sample);

                    [t, bf, ~, sd, pes, p_value] = stat.ttest_and_bf(first_sample, second_sample);
                    stat_data.time(sample, :)    = round(axis(sample));
                    stat_data.tValues(sample, :) = t;
                    stat_data.pValues(sample, :) = p_value;
                    stat_data.BFs(sample, :)     = bf;
                    stat_data.pes(sample, :)     = pes;
                    stat_data.sd(sample, :)      = sd;

                end
                path2save = [folder, filesep, 'EBR',  filesep, 'stat']; 
                if ~exist(path2save, 'file')
                    mkdir(path2save);
                end

                writetable(struct2table(stat_data), [path2save, filesep, contrasts{contranst_id}, '.csv']);
                data2plot = stat_data.BFs;
                data2plot(data2plot<3) = nan;

                cond1_id = find(strcmp(comp_names, ['c_' strrep(strrep(contrast_arr{1}, '&', 'x'), ' ', '_')]));
                cond2_id = find(strcmp(comp_names, ['c_' strrep(strrep(contrast_arr{2}, '&', 'x'), ' ', '_')]));

                colors = [0, 255, 0;
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

                line_pos =  min_val-abs(range*.1*contranst_id);
                line_dist = abs(range*.025);
                data2plot(~isnan(data2plot)) = line_pos;

                plot(axis, data2plot, 'Color', colors(cond1_id, :), 'LineWidth', 2,  'Marker','s', 'LineStyle','-',  'MarkerFaceColor', colors(cond1_id, :), 'MarkerSize', 1.5,  'HandleVisibility', 'off', 'Parent', fig);
                data2plot(~isnan(data2plot)) = line_pos-line_dist;
                plot(axis, data2plot, 'Color', colors(cond2_id, :), 'LineWidth', 2,  'Marker','s', 'LineStyle','-',  'MarkerFaceColor', colors(cond2_id, :), 'MarkerSize', 1.5,  'HandleVisibility', 'off', 'Parent', fig);
            
            end

            ylim([min_val-abs(range*.15*num_of_contrasts), max_val+0.1*range] );
            print_log(['Data has been analyzed successfully! ', num2str(toc) ' seconds'], log);    

        end

        
        function print_figure(blinks, condition_names, axis, fig)
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
            for condition = 1:length(condition_names)
% 
%                 data4max = blinks{condition};
%                 data4max = data4max(:, first_sample:first_sample+length(axis)-1);
%                 [val, pos] = max(data4max(:, 1001:3000)');
%                 max_data.pos(:, condition) = pos';
%                 max_data.val(:, condition) = val';

                [blink_avg, blink_se] = blinks_analysis.calculate_params(blinks{condition});
                blink_avg = blink_avg(1:length(axis));

                blink_se  = blink_se(1:length(axis));
                h = fill([axis fliplr(axis)], [blink_avg-blink_se  fliplr(blink_avg+blink_se)], cmap3(condition,:),  'LineStyle', ':', 'EdgeColor', cmap2(condition,:), 'HandleVisibility', 'off', 'Parent', fig);
                set(h,'facealpha',.5);

            end

            for condition = 1:length(condition_names)
                blink_avg = blinks_analysis.calculate_params(blinks{condition});
                blink_avg = blink_avg(1:length(axis));

                plot(axis, blink_avg, 'Color', cmap2(condition, :), 'DisplayName', condition_names{condition}, 'Parent', fig);
            end
            xtickformat('%,.4g');

   
            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('EBR', 'FontWeight','bold');

            xlim([axis(1), axis(end)]);
            legend('Location', 'Best');      

        end
        
        function [blink_avg, blink_se] = calculate_params(blinks)
            n  = size(blinks, 1);
            blink_avg  = mean(blinks);
            blink_se = std(blinks)./(n.^0.5)';
        end
        function exit(~, ~) % exit
            close(findobj('type', 'figure', 'name', 'CHAP - Eye-blinks Analysis'));
        end

        function data = show_statistical_vars(data, figure)
            hp = uipanel('Parent',figure, 'Units', 'Pixels', 'BackgroundColor', [0.1 0.32 0.46], 'BorderWidth', 0, 'Position', [33 320 300 350]);
            title_pos_x = 10;
            edit_pos_x  = 195;
            edit_pos_y  = 300;

            gui_lib.uicontrol_text(hp, 'To:', [title_pos_x edit_pos_y-40 50 30]);
            tooltip = '<html><b>To</b><br><i>Select last event for analyze</html>';

            range = round(data.x_axis);

            to = uicontrol(hp,'Style','popupmenu',...
                'Position',[edit_pos_x-40 edit_pos_y-40 120 30], 'TooltipString', tooltip, 'String', range(2:end), 'Value', size(range, 2)-1);%samples pre-event


            gui_lib.uicontrol_text(hp, 'From:', [title_pos_x edit_pos_y 50 30]);
            tooltip = '<html><b>From</b><br><i>Select first event for analyze</html>';
            from = uicontrol(hp,'Style','popupmenu',...
                'Position',[edit_pos_x-40 edit_pos_y 120 30], 'TooltipString', tooltip, 'String', range(1:end-1), 'Value', 1, 'Callback',{@blinks_analysis.change_to, to, range});%samples pre-event



            gui_lib.uicontrol_text(hp, 'Parameter:', [title_pos_x edit_pos_y-80 90 30]);
            tooltip = '<html><b>Parameter</b><br><i>Select the desired parameter for analysis</html>';
            measure = uicontrol(hp,'Style','popupmenu',...
                'Position',[edit_pos_x-40 edit_pos_y-80 120 30], 'TooltipString', tooltip, 'String', [{'mean'}; {'peak'};{'peak latency'}; {'dip'};{'dip latency'};], 'Value', 1);

            
            data.configuration.var.from        = from;
            data.configuration.var.to          = to;
            data.configuration.var.measure     = measure;

        end
        function change_to(src, ~, to, events)
            val_from = gui_lib.get_popupmenu_val(src);
            val_to   = gui_lib.get_popupmenu_val(to);

            if(isempty(str2num(char(val_from))))
                id       = find((strcmp(events, val_from))); 
            else
                id       = find(events==str2num(val_from)); 
            end
            events   = events(id+1:end);
            set(to, 'String', events);
            if(isempty(str2num(char(val_to))))
                to_id = find((strcmp(events, val_to)));
            else
                to_id = find(events==str2num(val_to)); 
            end
            if(to_id)
                set(to, 'Value', to_id);
            else
                set(to, 'Value', size(events,1));
            end
        end
        
        function stat_data = calc_stat_params(~, ~, data, blinks, figure)
            val_from = str2num(char(gui_lib.get_popupmenu_val(data.configuration.var.from)));
            from_id = find(round(data.x_axis)==val_from);
            val_to = str2num(char(gui_lib.get_popupmenu_val(data.configuration.var.to)));
            to_id = find(round(data.x_axis)==val_to);
            
            measure = char(gui_lib.get_popupmenu_val(data.configuration.var.measure));

            comd_names = data.configuration.comp_names;
            for condition = 1:length(comd_names)
                condition_data = blinks{condition};
                avgs.(char(comd_names(condition))) = mean(condition_data(:, from_id:to_id), 2, 'omitnan');
                if strcmp(measure,'peak') || strcmp(measure,'peak latency')
                    [val, pos]= max(condition_data(:, from_id:to_id)');

                    if strcmp(measure,'peak')
                        avgs.(char(comd_names(condition))) = val';
                    else
                        avgs.(char(comd_names(condition))) = round(data.x_axis(from_id+pos-1))';
                    end
                elseif strcmp(measure,'dip') || strcmp(measure,'dip latency')
                    [val, pos]= min(condition_data(:, from_id:to_id)');

                    if strcmp(measure,'dip')
                        avgs.(char(comd_names(condition))) = val';
                    else
                        avgs.(char(comd_names(condition))) = round(data.x_axis(from_id+pos-1))';
                    end
                end
            end

            conp_names_fixed = cellfun(@(x) x(3:end), comd_names, 'UniformOutput', false);
            conp_names_fixed = strrep(strrep(conp_names_fixed, '_x_', ' & '),'_',' ');

            
            % descriptive
            avgs_mat            = struct2array(avgs);
            num_of_subjects     = size(avgs_mat, 1);
            num_of_conditions   = size(avgs_mat, 2);
            means_conditions    = mean(avgs_mat, 1, 'omitnan');
            stds_conditions     = std(avgs_mat, 'omitnan');
            stes_conditions     = stds_conditions/sqrt(num_of_subjects);
            ts                  = tinv(0.975, num_of_subjects-1)*stes_conditions;
            CI_lower            = means_conditions-ts;
            CI_upper            = means_conditions+ts;
            
            width = 822;
            hight = 150;
            hpt   = uipanel('Parent',figure, 'Units','Pixels', 'BackgroundColor',[0.1 0.32 0.46], 'Position',[38 90 width  150]);
            y     = 147-hight;
            x     = 1;
            stat_data.output.measure = measure;
            stat_data.output.data    = struct2table(avgs);

            stat_data.output.descriptive = table(conp_names_fixed, (num_of_subjects.*ones(1, num_of_conditions))', means_conditions', stds_conditions', stes_conditions', CI_lower', CI_upper', 'VariableNames',{'Condition' 'N', 'Mean', 'SD', 'SE', 'CI_lower', 'CI_upper'});
            data_table = table2cell(stat_data.output.descriptive);
            header2 = stat_data.output.descriptive.Properties.VariableNames;
            uitable(hpt, 'Data', data_table, 'position', [x y width-1 hight-1],...
                      'ColumnName', header2', 'ColumnWidth',{370, 60, 70, 70, 70, 70, 70}); 
            guidata(figure, stat_data)
            
        end
        function save_stat(figure, folder)
            data = guidata(figure);
            output = data.output;
            project_name  = {''}; %prepare string
            prompt1         = {'Please enter analyzing name'};
            dlg_title1      = 'Please enter analyzing name'; 
            num_lines       = 1; %Num of lines
            answer1         = inputdlg(prompt1, dlg_title1, num_lines, project_name); %Create window
            
            if sum(size(answer1))>0 && ~strcmp(answer1{1}, '')
                project_name = answer1{1};
            else
                return;
            end

            writetable(output.descriptive, strcat(folder, filesep, 'EBR', filesep, project_name, '_', output.measure, '.csv'));
            writetable(output.data, strcat(folder, filesep, 'EBR', filesep, project_name, '_', output.measure, '_full.csv'));
        end
    end
end
