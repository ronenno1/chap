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
            fig = figure('Name', 'idttest');
            hold off
            plot(axis, stat_data.BFs, 'DisplayName', '');
            above_3 = stat_data.BFs;
            above_3(above_3<3) = NaN;
            hold on
            plot(axis, above_3,  'Color', [1, 0, 0], 'DisplayName', '');
            xtickformat('%,.4g');
            title(condition_str);
            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('BF_{10}', 'FontWeight','bold');
            xlim([-1000, 9000]);
            set(gca,'FontWeight','bold');

            max_val = max(stat_data.BFs);
            min_val = min(stat_data.BFs);

            range = max_val-min_val;

            ylim([0, max_val+0.1*range] )

            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs10.png']);
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs10.fig']);


            hold off
            %% BF01
            plot(axis, 1./stat_data.BFs, 'DisplayName', '');

            above_3 = 1./stat_data.BFs;
            above_3(above_3<3) = NaN;
            hold on

            plot(axis, above_3,  'Color', [1, 0, 0], 'DisplayName', '');

            xtickformat('%,.4g');
            title(condition_str);
            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('BF_{01}', 'FontWeight','bold');
            xlim([-1000, 9000]);

            set(gca,'FontWeight','bold');

            max_val = max(1./stat_data.BFs);
            min_val = min(1./stat_data.BFs);

            range = max_val-min_val;

            ylim([0, max_val+0.1*range] )

            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs01.png']);
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '_bfs01.fig']);

            hold off
            %% All data
            avg1 = nanmean(full_data1);
            avg2 = nanmean(full_data2);

            x = axis(1:min_axis)';
            d = avg1(1:min_axis)';
            scattering_data1 = nanstd(full_data1)./(size(full_data1,1).^0.5);
            scattering_data1 = scattering_data1(1:min_axis)';
            h = fill([x;flipud(x)],[d-scattering_data1;flipud(d+scattering_data1)], [1, 0, 0], 'LineStyle', '-', 'EdgeColor', [1, 0, 0], 'HandleVisibility', 'off');
            set(h,'facealpha',.2)

            hold on
            d = avg2(1:min_axis)';
            scattering_data2 = nanstd(full_data2)./(size(full_data2,1).^0.5);
            scattering_data2 = scattering_data2(1:min_axis)';
            h = fill([x;flipud(x)],[d-scattering_data2;flipud(d+scattering_data2)], [0, 1, 0], 'LineStyle', '-', 'EdgeColor', [0, 1, 0], 'HandleVisibility','off');
            set(h,'facealpha',.2)

            a = plot(axis(1:min_axis), avg1(1:min_axis),  'Color', [1, 0, 0], 'DisplayName', first_name);

            b = plot(axis(1:min_axis), avg2(1:min_axis),  'Color', [0, 1, 0], 'DisplayName', second_name);

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

            plot(axis(1:min_axis), BF10(1:min_axis), 'Color', [.75, .75, .75], 'LineWidth', 2, 'DisplayName', '');


            BF01 = 1./stat_data.BFs;
            BF01(BF01<3) = NaN;
            BF01(~isnan(BF01)) = line_pos;

            plot(axis(1:min_axis), BF01(1:min_axis), 'Color', [0, 0, 0], 'LineWidth', 2, 'DisplayName', '');

            % legend([a, b], {'Task', 'Ctrl'}, 'Location', 'Best');

            xlabel('Time [ms]', 'FontWeight','bold');
            ylabel('Relative Pupil Size [%]', 'FontWeight','bold');
            % ylabel('\Delta Relative Pupil Size [%]', 'FontWeight','bold');

            xtickformat('%,.4g');
            title(condition_str);
            xlim([-1000, 9000]);
            set(gca,'FontWeight','bold');

            legend('Location', 'Best');



            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '.png']);
            output.save_figure2(get(fig,'Children'), [output_folder_name, filesep, condition_str '.fig']);


            close(fig);
            save(['stat_data_' condition_str '.mat'], 'stat_data');
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
            for i=required_condition
                min_length = min(min_length, size(group_all_data.(full_condition_names{required_condition(i)}).data, 2));
            end
            for i=required_condition
                total = total + group_all_data.(full_condition_names{required_condition(i)}).data(:, 1:min_length);
            end
            data = total./length(required_condition);
            axis = group_all_data.x_axis;    
        end
        
        
        % function [t, bf, N, sd, pes, p_value] = ttest_and_bf(data1, data2, independent)
%     [t, sd, N, v] = my_ttest(data1, data2, independent);
%     if t == Inf || isnan(t)
%         t       = NaN;
%         bf      = NaN;
%         pes     = NaN;
%         p_value = NaN;
%         return;
%     end
%     r   = 0.707;
%     a   = (1+(t^2)/v)^-((v+1)/2);
%     fun = @(g) ((1+N*g*(r^2)).^(-1/2)) .* ((1+(t.^2)./((1+N.*g*(r^2)).*v)).^-((v+1)./2)) .* ((2.*pi).^(-1/2)) .* (g.^(-3/2)) .* exp(-1./(2.*g));
%     b   = integral(fun, 0, inf);
%     
%     bf      = b/a;
%     pes     = t^2/(t^2 + v);
%     p_value = 1-fcdf(t^2, 1, v);
% end

% function [t, sd, N, v] = my_ttest(data1, data2 , independent)
% 
%    if ~independent
%        data = data1- data2; 
%        d  = nanmean(data);
%        sd = nanstd(data);
%        N = sum(~isnan(data), 1);
%        se = sd/(N^0.5);
%        t = d/se;
%        v = N-1;
%         return;
%    end
%    sd = -1;
%    d1 = nanmean(data1);
%    d2 = nanmean(data2);
%    ss = sum((d1-data1).^2) + sum((d2-data2).^2);
%    n1 = sum(~isnan(data1), 1);
%    n2 = sum(~isnan(data2), 1);
%    N  = n1*n2/(n1+n2);
%    v  = n1 + n2 - 2;
%    s2 = ss/v;
%    t  = (d1-d2)/(sqrt((s2/n1) + (s2/n2)));
% end



    end
end