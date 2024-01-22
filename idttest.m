classdef idttest
    methods(Static)

        function run_idttest()
            debug = 0;
            if ~debug
                comparison_info  = {'', ''}; 

                prompt1         = {'First condition', 'Second condition'};
                dlg_title1      = 'Condition names ';
                num_lines       = 1;
                answer1         = inputdlg(prompt1, dlg_title1, num_lines, comparison_info);
                if sum(size(answer1))>0 && ~strcmp(answer1{1}, '') && ~strcmp(answer1{2}, '')
                    first_name  = answer1{1};
                    second_name = answer1{2};
                else
                    disp('Error: Condition names have to be chosen');
                    return;
                end


                output_folder_name = uigetdir('select destination directory');
                if(~output_folder_name)
                    return;
                end


                [full_data1, axis1] = idttest.load_file();
                if ~axis1
                    return
                end

                [full_data2, axis2] = idttest.load_file();
            end
            condition_str = [first_name , ' vs. ' second_name];


            min_axis = min([size(full_data1, 2), size(full_data2, 2), length(axis1), length(axis2)]);

            full_data1 = full_data1(:, 1:min_axis);
            full_data2 = full_data2(:, 1:min_axis);

            axis = axis1(1: min_axis);


            stat_data = [];
            for sample = 1:min_axis
                data1 = full_data1(:, sample);
                data2 = full_data2(:, sample);
                [t, bf, n, sd, pes, pvalue]        = stat.ttest_and_bf(data1, data2, true);
                stat_data.se(:, sample)  = sd./(n.^0.5);
                stat_data.pes(:, sample) = pes;
                stat_data.tValues(:, sample)  = t;
                stat_data.pValues(:, sample)  = pvalue;
                stat_data.BFs(:, sample)      = bf;
            end

            %% BF10
            fig      = figure('Name', 'idttest');
            fig_axes = axes('Parent',fig);

            hold off
            plot(axis, stat_data.BFs, 'DisplayName', '', 'Parent', fig_axes);
            above_3 = stat_data.BFs;
            above_3(above_3<3) = NaN;
            above_3_g = stat_data.BFs;
            between_3_g = stat_data.BFs;
%             between_3_g(between_3_g>=3 | between_3_g<=1/3) = NaN; % It is not necessary to remove these values since it will be presented in a different color

            above_3_g(above_3_g<3 & above_3_g>1/3) = NaN;
            
            hold on
            plot(axis, above_3,  'Color', [1, 0, 0], 'DisplayName', '', 'Parent', fig_axes);
            xtickformat('%,.4g');
            title(condition_str);
            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('BF_{10}', 'FontWeight','bold');
            xlim([axis(1), axis(end)]);
            set(gca,'FontWeight','bold');

            max_val = max(stat_data.BFs);
            min_val = min(stat_data.BFs);

            range = max_val-min_val;

            ylim([0, max_val+0.1*range] )
            
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs10.png']);
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs10.fig']);


            hold off
            %% BF01
            plot(axis, 1./stat_data.BFs, 'DisplayName', '', 'Parent', fig_axes);

            above_3 = 1./stat_data.BFs;
            above_3(above_3<3) = NaN;

            hold on

            plot(axis, above_3,  'Color', [1, 0, 0], 'DisplayName', '', 'Parent', fig_axes);

            xtickformat('%,.4g');
            title(condition_str);
            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('BF_{01}', 'FontWeight','bold');
            xlim([axis(1), axis(end)]);

            set(gca,'FontWeight','bold');

            max_val = max(1./stat_data.BFs);
            min_val = min(1./stat_data.BFs);

            range = max_val-min_val;

            ylim([0, max_val+0.1*range] )

            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs01.png']);
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs01.fig']);

            hold off
            %% All data
            avg1 = mean(full_data1, 'omitnan');
            avg2 = mean(full_data2, 'omitnan');

            x = axis(1:min_axis)';
            d = avg1(1:min_axis)';
            scattering_data1 = std(full_data1, 'omitnan')./(size(full_data1,1).^0.5);
            scattering_data1 = scattering_data1(1:min_axis)';
            h = fill([x;flipud(x)],[d-scattering_data1;flipud(d+scattering_data1)], [1, 0, 0], 'LineStyle', '-', 'EdgeColor', [1, 0, 0], 'HandleVisibility', 'off');
            set(h,'facealpha',.2)

            hold on
            d = avg2(1:min_axis)';
            scattering_data2 = std(full_data2, 'omitnan')./(size(full_data2,1).^0.5);
            scattering_data2 = scattering_data2(1:min_axis)';
            h = fill([x;flipud(x)],[d-scattering_data2;flipud(d+scattering_data2)], [0, 1, 0], 'LineStyle', '-', 'EdgeColor', [0, 1, 0], 'HandleVisibility','off');
            set(h,'facealpha',.2)

            a = plot(axis(1:min_axis), avg1(1:min_axis),  'Color', [1, 0, 0], 'DisplayName', first_name, 'Parent', fig_axes);

            b = plot(axis(1:min_axis), avg2(1:min_axis),  'Color', [0, 1, 0], 'DisplayName', second_name, 'Parent', fig_axes);

            %% stat

            fig_data = get(fig, 'children');
            min_val = fig_data.YAxis.Limits(1);
            max_val = fig_data.YAxis.Limits(2);

            range = max_val-min_val;
            line_pos =  min_val-abs(range*.05);

            ylim([min_val-abs(range*.2), max_val+0.05*range] )

            BF10 = stat_data.BFs;
            BF10(BF10<3) = NaN;
            BF10(~isnan(BF10)) = line_pos;

            plot(axis(1:min_axis), BF10(1:min_axis), 'Color', [.75, .75, .75], 'LineWidth', 2, 'DisplayName', '', 'Parent', fig_axes);


            BF01 = 1./stat_data.BFs;
            BF01(BF01<3) = NaN;
            BF01(~isnan(BF01)) = line_pos;

            plot(axis(1:min_axis), BF01(1:min_axis), 'Color', [0, 0, 0], 'LineWidth', 2, 'DisplayName', '', 'Parent', fig_axes);

            % legend([a, b], {'Task', 'Ctrl'}, 'Location', 'Best');

            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('Relative Pupil Size [%]', 'FontWeight','bold');
            % ylabel('\Delta Relative Pupil Size [%]', 'FontWeight','bold');

            xtickformat('%,.4g');
            title(condition_str);
            xlim([x(1), x(end)]);
            set(gca,'FontWeight','bold');

            legend('Location', 'Best');



            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '.png']);
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '.fig']);


            close(fig);

            fig = figure('Name', 'idttest');
            fig_axes = axes('Parent',fig);

            plot(axis(1:min_axis), log10(between_3_g(1:min_axis)), 'Color', [0, 0, 1], 'LineWidth', 2, 'DisplayName', '', 'Parent', fig_axes);
            hold on
            plot(axis(1:min_axis), log10(above_3_g(1:min_axis)),  'Color', [1, 0, 0], 'LineWidth', 2, 'DisplayName', '', 'Parent', fig_axes);

          
            fig_data = get(fig, 'children');
            min_val = fig_data.YAxis.Limits(1);
            max_val = fig_data.YAxis.Limits(2);

            range = max_val-min_val;
            
            ylim([min_val-abs(range*.05), max_val+0.05*range] )
            lables = round(10*(10.^fig_data.YTick))/10;
            new_lables = arrayfun(@(x) stat.log_round(x), lables);
            new_lables = unique(new_lables);

            fig_data.YTick = log10(new_lables);
            fig_data.YTickLabel = new_lables;

            
            high_threshold = ones(min_axis, 1) * log10(3);

            plot(axis(1:min_axis), high_threshold, '--' , 'color', [0, 0, 0],  'LineWidth', 2, 'Parent', fig_axes);

            low_threshold = ones(min_axis, 1) * log10(1/3);

            plot(axis(1:min_axis), low_threshold, '--' , 'color', [0, 0, 0],  'LineWidth', 2, 'Parent', fig_axes);
            
            
            plot(axis(1:min_axis), ones(min_axis, 1)*log10(1), '--' , 'color', [0.75, 0.75, 0.75],  'LineWidth', 2, 'Parent', fig_axes);


            
            xtickformat('%,.4g');
            title(condition_str);
            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('BF_{10}', 'FontWeight','bold');
            xlim([axis(1), axis(end)]);

            set(gca,'FontWeight','bold');

            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs.png']);
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs.fig']);
            
            save([output_folder_name, filesep 'stat_data_' condition_str '.mat'], 'stat_data');
        end
        
        
        function [data, axis] = load_file()
            data = 0;
            axis = 0;
            [file_name, path]  = uigetfile('*.mat', 'Select mat file ');
            if ~path
               disp('Error: *.mat file have to be chosen');
               return;
            end

            full_mat_name = [path, file_name];

            group = load(full_mat_name);
            group_all_data = group.total_data;
            full_condition_names = group_all_data.configuration.comp_names;
            condition_names = cellfun(@(x) x(3:end), full_condition_names, 'UniformOutput', false);
            condition_names = strrep(strrep(condition_names, '_x_', ' & '), '_', ' ');
            required_condition        = listdlg('PromptString',...
                                       'Select condition:',...
                                       'ListSize',[200, 100],...
                                       'SelectionMode','multiple',...
                                       'ListString', condition_names);    

            if isempty(required_condition)
               disp('Error: Conditions have to be chosen');
               return;
            end
            
            total = 0;
            min_length = inf;
            if length(required_condition) >1
                for i=required_condition
                    min_length = min(min_length, size(group_all_data.(full_condition_names{i}).data, 2));
                end
            else
                min_length = size(group_all_data.(full_condition_names{required_condition}).data, 2);
            end
            
            if length(required_condition) >1

                for i=required_condition
                    total = total + group_all_data.(full_condition_names{i}).data(:, 1:min_length);
                end
            else
                total = group_all_data.(full_condition_names{required_condition}).data(:, 1:min_length);
            end

            data = total./length(required_condition);
            axis = group_all_data.x_axis;    
        end

    end
end