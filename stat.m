classdef stat
    methods(Static)

        function configuration = statistical_vals(configuration)
            configuration.vals.measure      = gui_lib.get_popupmenu_val(configuration.var.measure);
            configuration.vals.approach     = gui_lib.get_popupmenu_val(configuration.var.approach);
            configuration.vals.output_type  = gui_lib.get_popupmenu_val(configuration.var.output_type);
            configuration.vals.from         = str2num(gui_lib.get_popupmenu_val(configuration.var.from));
            configuration.vals.to           = str2num(gui_lib.get_popupmenu_val(configuration.var.to));
        end
        
        function calc_and_plot(src, evnt, fig, figure, total_data, statistical_data, comp_names, next, back, configuration)

            configuration = stat.statistical_vals(configuration);

            statistical_data = stat.do_stat(total_data, comp_names, configuration, statistical_data);
            statistical_data.configuration = configuration;
            stat.plot_statistics(0,0, fig, figure, statistical_data, total_data, 1, 0, next, back);
        end

        function plot_statistics(src, evnt, fig, figure, statistical_data, total_data, contrast_id, action, next, back)
            statistical_data.configuration = stat.statistical_vals(statistical_data.configuration);
            bayesian = strcmp(statistical_data.configuration.vals.approach, 'Bayesian');
            descriptive = strcmp(statistical_data.configuration.vals.output_type, 'descriptive');

            approach_str = 'BFs';

            if ~bayesian
                approach_str = 'pValues';
            end
            contrast_names = fieldnames(statistical_data.contrasts.(approach_str));
            header = statistical_data.contrast_names;
            
            if (contrast_id+action)>size(contrast_names, 1) || contrast_id+action<1
                return;
            end
            
            bayesian_contrast_data       = statistical_data.contrasts.BFs.(char(contrast_names(contrast_id+action)));
            classical_contrast_data_bonf = statistical_data.contrasts.pValues_bonf.(char(contrast_names(contrast_id+action)));
            contrast_data                = bayesian_contrast_data;
            
            if ~bayesian
                contrast_data = classical_contrast_data_bonf;
            end
            
            
            contrast_data_effect  = statistical_data.contrasts.pes.(char(contrast_names(contrast_id+action)));
          
            x_axis = statistical_data.configuration.x_axis;
            if (size(statistical_data.configuration.x_axis, 2)> size(contrast_data, 2))
                x_axis = x_axis(1:size(contrast_data, 2));
            end
            data2plot = contrast_data(1:size(x_axis, 2));
            max_bf = max(data2plot);
            if bayesian
                data2plot(data2plot<3) = nan;
                data2plot(data2plot>=3) = -max_bf/50;
            else
                data2plot_001(data2plot>.001) = nan;
                data2plot_001(data2plot<.001) = -max_bf/50;
                data2plot_01(find(data2plot>.01 | data2plot<.001)) = nan;
                data2plot_01(find(data2plot<.01 & data2plot>=.001)) = -max_bf/50;
                
                data2plot_05(find(data2plot>.05 | data2plot<.01)) = nan;
                data2plot_05(find(data2plot<.05 & data2plot>=.01)) = -max_bf/50;
            end
            
            data2plot_ns = contrast_data(1:size(x_axis, 2));
            data2plot_ns(data2plot_ns>(1/3)) = nan;
            data2plot_ns(data2plot_ns<=(1/3)) = -max_bf/50;

            cla(fig);
            plot(fig, x_axis, contrast_data(1:size(x_axis, 2)), 'blue', 'LineWidth',2);

            hold on
            
            if bayesian
                plot(fig, x_axis, 1./contrast_data(1:size(x_axis, 2)), 'red', 'LineWidth',2);
%                 grouper.do_plot(total_data, total_data.configuration, total_data.configuration.rate, fig);
                
                plot(fig, x_axis, data2plot(1:size(x_axis, 2)), 'black', 'LineWidth', 0.1, 'Marker','s', 'LineStyle','none',  'MarkerFaceColor', 'black');
                plot(fig, x_axis, data2plot_ns(1:size(x_axis, 2)), 'green', 'LineWidth', 0.1, 'Marker','s', 'LineStyle','none',  'MarkerFaceColor', 'green');
                ylim([-inf inf])
            else
                plot(fig, x_axis, data2plot_001(1:size(x_axis, 2)), 'color', [0, 0, 0], 'LineWidth', 0.1, 'Marker','s', 'LineStyle','none',  'MarkerFaceColor', [0, 0, 0]);
                plot(fig, x_axis, data2plot_01(1:size(x_axis, 2)), 'color', [.3, .3, .3], 'LineWidth', 0.1, 'Marker','s', 'LineStyle','none',  'MarkerFaceColor', [.3, .3, .3]);
                plot(fig, x_axis, data2plot_05(1:size(x_axis, 2)), 'color', [.6, .6, .6], 'LineWidth', 0.1, 'Marker','s', 'LineStyle','none',  'MarkerFaceColor', [.6, .6, .6]);

            end
            set(fig,'FontWeight','bold');

            xlabel(fig, 'Time [ms]', 'FontWeight','bold');
            

            if bayesian
                ylabel(fig, 'BF', 'FontWeight','bold');
                legend(fig, ['BF_{10}'; 'BF_{01}']);
            else
                ylabel(fig, 'p-value', 'FontWeight','bold');
                legend(fig, 'p-value');

            end
            xlim(fig, [statistical_data.configuration.x_axis(1) statistical_data.configuration.x_axis(end)])

            title(fig, header(contrast_id+action));        
            set(next, 'Visible', 'off');
            set(back, 'Visible', 'off');
            if (contrast_id+action)<size(contrast_names, 1)
                set(next, 'Visible', 'on');
            end
            if contrast_id+action>1
                set(back, 'Visible', 'on');
            end
            set(next, 'callback', {@stat.plot_statistics, fig, figure, statistical_data, total_data, contrast_id+action, 1, next, back}); 
            set(back, 'callback', {@stat.plot_statistics, fig, figure, statistical_data, total_data, contrast_id+action, -1, next, back}); 
            width = 822;
            hight = 150;
    
            hpt   = uipanel('Units','Pixels', 'BackgroundColor',[0.1 0.32 0.46], 'Position',[38 90 width  150]);
            y     = 147-hight;
            x     = 1;


            header2 = [cellstr('Comparison'), cellstr('t'), cellstr('BF 10'), cellstr('BF 01'), cellstr('Cohen''s d'), 'N'];
            data_table = cell(size(header, 1), 6);
            
            for i=1:size(header, 1)
                comparison = header(i);
                BF_10 = statistical_data.contrasts.total.bf(i);
                cohensd = statistical_data.contrasts.total.cohensd(i);

                t = statistical_data.contrasts.total.t(i);
                N = statistical_data.contrasts.total.N(i);
                data_table(i,:) = [cellstr(comparison), t, BF_10, 1./BF_10, cohensd, N];
            end
            
            uitable(hpt, 'Data', data_table, 'position', [x y width-1 hight-1],...
                          'ColumnName', header2', 'ColumnWidth',{500,70, 50, 50, 70, 40}); 

            if ~bayesian
                data_table = table2cell(statistical_data.output.clasical);
                header2 = statistical_data.output.clasical.Properties.VariableNames;
                uitable(hpt, 'Data', data_table, 'position', [x y width-1 hight-1],...
                          'ColumnName', header2', 'ColumnWidth',{500,35, 35, 50, 50, 70, 40}); 
        

            end
            if descriptive
                data_table = table2cell(statistical_data.output.descriptive);
                header2 = statistical_data.output.descriptive.Properties.VariableNames;
                uitable(hpt, 'Data', data_table, 'position', [x y width-1 hight-1],...
                          'ColumnName', header2', 'ColumnWidth',{370, 60, 70, 70, 70, 70, 70}); 
            end
            tooltip = '<html><i>Save data';
            gui_lib.uicontrol_button(figure, [162 40 155 30],'Save data', {@(~,~)stat.save_stat(statistical_data.output, statistical_data.configuration.paths)}, tooltip); 
                     
        end
        
        function save_stat(output, paths)
            folder_name  = {''}; %prepare string
            prompt1         = {'Please enter analyzing name'};
            dlg_title1      = 'Please enter analyzing name'; 
            num_lines       = 1; %Num of lines
            answer1         = inputdlg(prompt1, dlg_title1, num_lines, folder_name); %Create window
            
            if sum(size(answer1))>0 && ~strcmp(answer1{1}, '')
                folder_name = answer1{1};
            else
                return;
            end
            path2save = strcat(paths.stat_output_folder_name, filesep, folder_name);
             if ~exist(path2save, 'dir')
                mkdir(path2save);
            end

            writetable(output.descriptive, strcat(path2save, filesep, 'descriptive.csv'));
            writetable(output.inference, strcat(path2save, filesep, 'inference.csv'));
            
            writetable(output.bayesian_sequential, strcat(path2save, filesep, 'bayesian_sequential.csv'));
            writetable(output.raw_data, strcat(path2save, filesep, 'avgs_table.csv'));

        end
        
        function [t, sd, N, v] = my_ttest(data1, data2)
           data = data1- data2; 
           d  = nanmean(data);
           sd = nanstd(data);
           N = sum(~isnan(data), 1);
           se = sd/(N^0.5);
           t = d/se;
           v = N-1;         
        end
        
        function [t, bf, N, sd, pes, p_value] = ttest_and_bf(data1, data2)
            [t, sd, N, v] = stat.my_ttest(data1, data2);
            if t == Inf || isnan(t)
                t = NaN;
                bf = NaN;
                pes = NaN;
                p_value = NaN;
                return;
            end
            r   = 0.707;
            a   = (1+(t^2)/v)^-((v+1)/2);
            fun = @(g) ((1+N*g*(r^2)).^(-1/2)) .* ((1+(t.^2)./((1+N.*g*(r^2)).*v)).^-((v+1)./2)) .* ((2.*pi).^(-1/2)) .* (g.^(-3/2)) .* exp(-1./(2.*g));
            b   = integral(fun, 0, inf);
            
            bf = b/a;
            pes = t^2/(t^2 + v);
            p_value = 1-fcdf(t^2, 1, v);
            
        end
        
        function num = rounder(X, N)
        if (~exist('N', 'var'))
                N = 3;
            end
            num = round(X*(10^N))/10^N;
        end

        function [sse, ssa, sst, sss] = calc_ss(avgs_table_mat)
            num_of_subjects   = size(avgs_table_mat, 1);
            num_of_conditions = size(avgs_table_mat, 2);
            means_conditions = nanmean(avgs_table_mat);
            grand_total      = nanmean(means_conditions);
            extended_means_table = repmat(grand_total, 1, num_of_conditions);
            diff_conds           = means_conditions-extended_means_table;
            diff_squere          = diff_conds.*diff_conds;
            
            ssa = sum(sum(diff_squere))*size(avgs_table_mat, 1);
            means_subjects = nanmean(avgs_table_mat, 2);
            extended_means_table = repmat(grand_total, num_of_subjects, 1);

            diff_subs = means_subjects-extended_means_table;
            diff_squere = diff_subs.*diff_subs;
            sss = num_of_conditions*sum(sum(diff_squere));
            
            %% sse
            extended_means_table = repmat(grand_total, size(avgs_table_mat));
            diff_total = avgs_table_mat-extended_means_table;
            diff_squere = diff_total.*diff_total;
            sst = sum(sum(diff_squere));
            sse = sst-sss-ssa;
        end
        
        function stat_data = do_stat(total_data, comp_names, configuration, statistical_data)
            stat_data.configuration.paths  = total_data.paths;            
            stat_data.configuration.x_axis = total_data.x_axis;
            if (~exist('statistical_data', 'var'))
                statistical_data = [];
            end
            if (~exist('configuration', 'var'))
                from = total_data.x_axis(1);
                to = total_data.x_axis(end);
                measure = 'mean';
                approach = 'Bayesian';
                output_type = 'descriptive';
            else
                from = configuration.vals.from;
                to = configuration.vals.to;
                measure = configuration.vals.measure;
                approach = configuration.vals.approach;
                output_type = configuration.vals.output_type;
            end
            bayesian = strcmp(approach, 'Bayesian');

            min_length = inf;
            for comp = 1:size(comp_names, 1)
                min_length = min(min_length, size(total_data.(char(comp_names(comp))).data, 2));
            end

            contrasts = nchoosek(1:size(comp_names, 1), 2);
            
            from_ana = 1; 
            if from > total_data.x_axis(1)
                from_ana = find(total_data.x_axis<=from, 1, 'last');
            end
            to_ana = length(total_data.x_axis)-1; 
            if to < total_data.x_axis(end)
                to_ana = find(total_data.x_axis>=to, 1, 'first')-1;
            end
            

            for comp = 1:size(comp_names, 1)
                full_data = total_data.(char(comp_names(comp))).data;
                avgs.(char(comp_names(comp))) = nanmean(full_data(:, from_ana:to_ana), 2);
                if strcmp(measure,'peak') || strcmp(measure,'peak latency')
                    [val, pos]= nanmax(full_data(:, from_ana:to_ana)');

                    if strcmp(measure,'peak')
                        avgs.(char(comp_names(comp))) = val';
                    else
                        avgs.(char(comp_names(comp))) = configuration.x_axis(from_ana+pos)';
                    end
                elseif strcmp(measure,'dip') || strcmp(measure,'dip latency')
                    [val, pos]= nanmin(full_data(:, from_ana:to_ana)');

                    if strcmp(measure,'dip')
                        avgs.(char(comp_names(comp))) = val';
                    else
                        avgs.(char(comp_names(comp))) = configuration.x_axis(from_ana+pos)';
                    end
                end
            end

            avgs_table = struct2table(avgs);
            avgs_table_mat = avgs_table{:,1:end};

            num_of_subjects   = size(avgs_table_mat, 1);
            num_of_conditions = size(avgs_table_mat, 2);
            if(num_of_subjects<2 || num_of_conditions<2)
                return;
            end

            comp_names_fixed = cellfun(@(x) x(3:end), comp_names, 'UniformOutput', false);
            comp_names_fixed = strrep(strrep(comp_names_fixed, '_x_', ' & '),'_',' ');

            
            % descriptive
            means_conditions = nanmean(avgs_table_mat, 1);
            stds_conditions  = nanstd(avgs_table_mat);
            stes_conditions  = stds_conditions/sqrt(num_of_subjects);
            ts = tinv(0.975, num_of_subjects-1)*stes_conditions;
            CI_lower = means_conditions-ts;
            CI_upper = means_conditions+ts;

            stat_data.output.descriptive = table(comp_names_fixed, (num_of_subjects.*ones(1, num_of_conditions))', means_conditions', stds_conditions', stes_conditions', CI_lower', CI_upper', 'VariableNames',{'Condition' 'N', 'Mean', 'SD', 'SE', 'CI_lower', 'CI_upper'});
            
            
            if ~ bayesian
                %% anova
                [sse, ssa] = stat.calc_ss(avgs_table_mat);
                dfa = num_of_conditions-1;
                dfe = (num_of_conditions-1)*(num_of_subjects-1);

                mse = sse / dfe;
                msa = ssa / dfa;

                F = msa/mse;
                p = 1-fcdf(F, dfa, dfe);
                pes = ssa/(ssa+sse);

                comps_t = -ones(num_of_conditions,num_of_conditions-1);
                for j=2:num_of_conditions
                    for i=2:num_of_conditions
                        if i==j
                           comps_t(i, j-1) = i-1;
                        elseif i>j
                            comps_t(i, j-1) = 0;
                        end
                    end
                end


                means_conditions_extended      = zeros(num_of_conditions-1, num_of_conditions);
                means_conditions_extended(1,:) = means_conditions;
                denominator  = 1./sum(comps_t.^2);
                numerator    = num_of_subjects*(means_conditions_extended * comps_t).^2;
                ss_comp_all  = numerator(1, :).*denominator;
                comp_sse = zeros(1, size(comps_t, 2));
                comp_mse = zeros(1, size(comps_t, 2));

                dfe_comp = (num_of_subjects-1);
                for comp = 1:size(comps_t, 2)
                    small_table = avgs_table_mat;
                    small_table (:, (comps_t(:, comp)==0)) = [];
                    comp_sse(comp) = stat.calc_ss(small_table);
                    comp_mse(comp) = stat.calc_ss(small_table)/dfe_comp;
                end

                F_comp_all  = (ss_comp_all./comp_mse);
                p_comp__all  = 1-fcdf(F_comp_all, 1, dfe_comp);

                fixed_contrast_names = cell(size(comps_t, 2), 1);
                for contrast = 1:size(comps_t, 2)
                    name = '';
                    for pos = 1:size(comps_t, 1)
                        weight = comps_t(pos, contrast);
                        if weight ==0
                            continue;
                        end
                        if sign(weight)>0
                            name = strcat(name, '+');
                        else
                            name = strcat(name, '-');
                        end
                        name = strcat(name, num2str(abs(weight)), '*', comp_names_fixed(pos));
                    end
                    fixed_contrast_names(contrast, :) = name;
                end

                p_str = stat.rounder(p);
                if(~p_str)
                    p_str = p;
                end
                
                stat_data.output.inference.Comperison = ' ';
                stat_data.output.inference.dfa = stat.rounder(dfa);
                stat_data.output.inference.dfe = stat.rounder(dfe);

                stat_data.output.inference.F       = stat.rounder(F);
                stat_data.output.inference.MSE     = mse;
                stat_data.output.inference.Pvalue  = p_str;
                stat_data.output.inference.pes     = {num2str(stat.rounder(pes))};

                stat_data.output.comps = table(fixed_contrast_names, ones(size(comps_t, 2), 1), dfe*ones(size(comps_t, 2), 1), F_comp_all', comp_mse', p_comp__all', cell(size(comps_t, 2), 1),...
                                               'VariableNames',{'Comperison', 'dfa', 'dfe', 'F', 'MSE', 'Pvalue', 'pes'});
                stat_data.output.clasical = [struct2table(stat_data.output.inference); stat_data.output.comps];  
            end

            stat_data.output.raw_data = avgs_table;

            stat_data.configuration.range = round(linspace(stat_data.configuration.x_axis(1), stat_data.configuration.x_axis(end), min_length));
            contrasts_table = cell(size(contrasts, 1));         
            percentages = zeros(1, 10);

            for contrast = 1:size(contrasts, 1)
                contrast_name = [char(comp_names(contrasts(contrast, 1))), '_vs_', char(comp_names(contrasts(contrast, 2)))];
                contrasts_table{contrast} = contrast_name;
                full_data1 = total_data.(char(comp_names(contrasts(contrast, 1)))).data;
                full_data2 = total_data.(char(comp_names(contrasts(contrast, 2)))).data;                

                total_avg1 = avgs.(char(comp_names(contrasts(contrast, 1))));
                total_avg2 = avgs.(char(comp_names(contrasts(contrast, 2))));
                
                [total_t, total_bf, N, total_sd, pes] = stat.ttest_and_bf(total_avg1, total_avg2);
                for s=1:size(total_avg1, 1)
                    [~, sequential_bf, ~, ~, ~] = stat.ttest_and_bf(total_avg1(1:s), total_avg2(1:s));
                    stat_data.contrasts.sequential.bf(contrast, s) = sequential_bf;
                end
                if isempty(statistical_data)
                    for sample = 1:min_length
                        percentage = round(100*((sample+min_length*(contrast-1))/(min_length*size(contrasts, 1))));
                        if percentage>0 && ~mod(percentage, 10) && percentages(percentage/10)==0 
                            percentages(percentage/10) = 1;
                        end

                        data1 = full_data1(:, sample);
                        data2 = full_data2(:, sample);
                        [t, bf, ~, ~, pes, pvalue] = stat.ttest_and_bf(data1, data2);
                        stat_data.contrasts.tValues.(char(['contrast_' num2str(contrast)]))(:, sample)  = t;
                        stat_data.contrasts.pValues.(char(['contrast_' num2str(contrast)]))(:, sample)  = pvalue;
                        stat_data.contrasts.BFs.(['contrast_' num2str(contrast)])(:, sample)      = bf;

                        mean_diff = nanmean(data1)-nanmean(data2); 
                        mean_std = nanmean(nanstd(data1)^2-nanstd(data2)^2); 
                        stat_data.contrasts.pes.(['contrast_' num2str(contrast)])(:, sample)      = pes;
                    end
                    pValues = stat_data.contrasts.pValues.(['contrast_' num2str(contrast)]);
                    num_of_samples           = length(pValues);
                    [sorted_data, sort_ids]  = sort(pValues);  
                    [~, real_ids]            = sort(sort_ids);  
                    factors                  = num_of_samples:-1:1;
                    sorted_pValues_bonf      = min(1, sorted_data.*factors);
                    sorted_pValues_bonf(2:num_of_samples) = max([sorted_pValues_bonf(1:num_of_samples-1); sorted_pValues_bonf(2:num_of_samples)]);
                    stat_data.contrasts.pValues_bonf.(['contrast_' num2str(contrast)])  = sorted_pValues_bonf(real_ids);
                   
                    num_of_samples         = length(total_data.x_axis);
                    data2save.bin          = (1:num_of_samples)';
                    data2save.x_axis       = total_data.x_axis';
                    data2save.tValues      = stat_data.contrasts.tValues.(['contrast_' num2str(contrast)])(1:num_of_samples)';
                    data2save.pValues      = stat_data.contrasts.pValues.(['contrast_' num2str(contrast)])(1:num_of_samples)';
                    data2save.pValues_bonf = stat_data.contrasts.pValues_bonf.(['contrast_' num2str(contrast)])(1:num_of_samples)';
                    data2save.BFs          = stat_data.contrasts.BFs.(['contrast_' num2str(contrast)])(1:num_of_samples)';
                    contrast_names_fixed   = strrep(strrep(strrep(contrast_name, 'c_', ''), '_x_', ' & '),'_',' ');
                    writetable(struct2table(data2save), [total_data.paths.stat_output_folder_name filesep contrast_names_fixed '.csv'])

                else
                    stat_data.contrasts.tValues      = statistical_data.contrasts.tValues;
                    stat_data.contrasts.pValues      = statistical_data.contrasts.pValues;
                    stat_data.contrasts.pValues_bonf = statistical_data.contrasts.pValues_bonf;
                    stat_data.contrasts.BFs          = statistical_data.contrasts.BFs;
                    stat_data.contrasts.pes          = statistical_data.contrasts.pes;
                end
                
                
                stat_data.contrasts.total.bf(contrast, :) = total_bf;
                stat_data.contrasts.total.t(contrast, :)  = total_t;
                stat_data.contrasts.total.N(contrast, :) = N;
                stat_data.contrasts.total.cohensd(contrast, :) = (nanmean(total_avg1)-nanmean(total_avg2))/total_sd;
            end
            all_fixed_contrast_names = cell(size(contrasts, 1), 1);
            for contrast = 1:size(contrasts, 1)
                all_fixed_contrast_names{contrast, :} = [char(comp_names_fixed(contrasts(contrast, 1))), ' - ', char(comp_names_fixed(contrasts(contrast, 2)))];
            end

            stat_data.output.bayesian_sequential = table(all_fixed_contrast_names,...
                                                         stat_data.contrasts.sequential.bf,...
                                                         'VariableNames', {'Comperison', 'BF10'});

          stat_data.output.inference = table(all_fixed_contrast_names,...
              stat_data.contrasts.total.t,...
              stat_data.contrasts.total.bf,...
              stat_data.contrasts.total.N,...
              stat_data.contrasts.total.cohensd,...
              'VariableNames', {'Comperison', 't', 'BF','N', 'cohensd'});
            stat_data.contrast_names = all_fixed_contrast_names;
        end 
    end
end