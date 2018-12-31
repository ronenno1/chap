classdef grouper    
    methods(Static)
        function [paths, files, overwrite] =  init_folders(file_type)
            paths     = [];
            files     = [];
            overwrite = 0;
            output_folder_name = uigetdir(['Select directory of ' file_type 's']);
            if(~output_folder_name)
                return;
            end
            
            project_name  = {''}; %prepare string
            prompt1         = {'Please enter project name'};
            dlg_title1      = 'Project name'; 
            num_lines       = 1; %Num of lines
            answer1         = inputdlg(prompt1, dlg_title1, num_lines, project_name); %Create window
            
            if sum(size(answer1))>0 && ~strcmp(answer1{1}, '')
                project_name = answer1{1}; %Log participant
            else
                return;
            end
            
            project_output_folder_name = output_folder_name;

            if sum(size(answer1))>0 && ~strcmp(answer1{1}, '')
                 project_name = answer1{1}; %Log participant
                 project_output_folder_name = strcat(output_folder_name, filesep, project_name);
                 if ~exist(project_output_folder_name, 'dir')
                    mkdir(project_output_folder_name);
                    overwrite = 1;
                 else
                     choice = questdlg('The project already exists. Do you want to overwrite it?', ...
                    'Overwrite', ...
                    'Yes', 'No', 'No');
                    % Handle response
                    switch choice
                        case 'Yes'
                            overwrite = 1;
                        case 'No'
                            overwrite = 0;
                    end
                 end
            end

            files = dir([output_folder_name filesep '*.' file_type]);
            bytes = {files.bytes}';
            bytes_arr = cell2mat(bytes);
            files = files(bytes_arr>10000);
            files = {files.name}';
            chap_output_folder_name = [output_folder_name filesep 'chap'];
            if ~exist(chap_output_folder_name, 'dir')
                mkdir(chap_output_folder_name);
            end
            
            png_output_folder_name = [project_output_folder_name filesep 'png'];
            if ~exist(png_output_folder_name, 'dir')
                mkdir(png_output_folder_name);
            end
            
            png_output_folder_name_err = [png_output_folder_name filesep 'outliers'];
            if ~exist(png_output_folder_name_err, 'dir')
                mkdir(png_output_folder_name_err);
            end
            
            
            stat_output_folder_name = [project_output_folder_name filesep 'stat'];
            if ~exist(stat_output_folder_name, 'dir')
                mkdir(stat_output_folder_name);
            end

            behave_output_folder_name = [project_output_folder_name filesep 'behave'];
            if ~exist(behave_output_folder_name, 'dir')
                mkdir(behave_output_folder_name);
            end

            
            fig_output_folder_name = [project_output_folder_name filesep 'fig'];
            if ~exist(fig_output_folder_name, 'dir')
                mkdir(fig_output_folder_name);
            end
            
            fig_output_folder_name_err = [fig_output_folder_name filesep 'outliers'];
            if ~exist(fig_output_folder_name_err, 'dir')
                mkdir(fig_output_folder_name_err);
            end
            
            
            mat_output_folder_name = [project_output_folder_name filesep 'mat'];
            if ~exist(mat_output_folder_name, 'dir')
                mkdir(mat_output_folder_name);
            end
            mat_output_folder_name_err = [project_output_folder_name filesep 'mat' filesep 'outliers'];
            if ~exist(mat_output_folder_name_err, 'dir')
                mkdir(mat_output_folder_name_err);
            end
            
            if(overwrite)
                delete([mat_output_folder_name filesep '*.mat']);
            end
            csv_output_folder_name = [project_output_folder_name filesep 'csv'];
            if ~exist(csv_output_folder_name, 'dir')
                mkdir(csv_output_folder_name);
            end
            paths.project_name               = project_name;
            paths.output_folder_name         = output_folder_name;
            paths.project_output_folder_name = project_output_folder_name;
            paths.chap_output_folder_name    = chap_output_folder_name;
            paths.png_output_folder_name     = png_output_folder_name;
            paths.png_output_folder_name_err = png_output_folder_name_err;

            paths.fig_output_folder_name     = fig_output_folder_name;
            paths.fig_output_folder_name_err = fig_output_folder_name_err;
            
            paths.behave_output_folder_name  = behave_output_folder_name;
            
            paths.mat_output_folder_name     = mat_output_folder_name;
            paths.mat_output_folder_name_err = mat_output_folder_name_err;

            paths.csv_output_folder_name     = csv_output_folder_name;
            
            paths.stat_output_folder_name    = stat_output_folder_name;
        end
        
        function total_data = do_avaraging(paths, comp_names)
            mat_files  = dir([paths.mat_output_folder_name filesep '*.mat']);
            bytes = {mat_files.bytes}';
            bytes_arr = cell2mat(bytes);
            mat_files = mat_files(bytes_arr>10000);

            if(isempty(mat_files))
                total_data = [];
                return;
            end
            mat_files  = {mat_files.name}';
            for file=1:size(mat_files, 1)
                mat_file = load([paths.mat_output_folder_name filesep char(mat_files(file))]);
                mat_file_data = mat_file.ploted_data;        
                for i = 1:size(comp_names, 1)
                    comp = char(comp_names(i));
                    if~isfield(mat_file_data, comp)
                        continue;
                    end
                    pupil_data = mat_file_data.(char(comp)).pupil';
                    if(file>1) %padding matrix
                        size_diff = size(total_data.(char(comp)).data, 2) - size(pupil_data, 2);
                        if(size_diff>0)
                            pupil_data = padarray(pupil_data, [0, size_diff], nan, 'post');
                        end
                        if(size_diff<0)
                            total_data.(char(comp)).data = padarray(total_data.(char(comp)).data, [0, -size_diff], nan, 'post');
                        end
                    end
                    total_data.(char(comp)).data(file, :) = pupil_data;
                    event_names = fieldnames(mat_file_data.(char(comp)).events);
                    for j=1:size(event_names, 1)
                        evnt = char(event_names(j));
                        total_data.(char(comp)).events.(evnt)(file,:) = mat_file_data.(char(comp)).events.(evnt);
                    end
                end
            end
        end
        
        function [total_data, comp_names_fixed] = do_plot(total_data, configuration, rate, fig)
            cla(fig);
            comp_names      = configuration.comp_names;
            PreEventNumber  = configuration.PreEventNumber_val;
            from           = configuration.from_val;
            Method         = configuration.Method_val;
            scattering     = configuration.scattering_val;
            relative       = configuration.relative_val;   

            ms   = round(1000/rate);
            
            hold off;
            Bins = ceil(max(configuration.BinsNumber_val, 1)/ms);
            
            for i = 1:size(comp_names, 1)
                comp = char(comp_names(i));
                lengths.(char(comp)) = floor(sum(sum(~isnan(total_data.(char(comp)).data)))/size(total_data.(char(comp)).data, 1));
            end
            
            cmap2 = [0, 255, 0;
                    255, 0, 0; 
                    0, 0, 255;
                    255, 0, 255;
                    0, 255, 255;
                    255, 255, 0; 
                    0, 0, 0]/255;
                
            cmap3 = cmap2+0.8;
            cmap3(cmap3>1) = 1;
            opengl_data = opengl('data');
                
            for comp_id = 1:size(comp_names, 1)
                comp = char(comp_names(comp_id));
                total_data.(char(comp)).avg = nanmean(total_data.(char(comp)).data, 1);
                total_data.(char(comp)).std = nanstd(total_data.(char(comp)).data, 1);

                event_names = fieldnames(total_data.(char(comp)).events);
                for j=1:size(event_names, 1)
                    evnt = char(event_names(j));
                    total_data.(char(comp)).avg_events.(evnt) = nanmean(total_data.(char(comp)).events.(evnt));
                end
                
                total_data.(char(comp)).avg = total_data.(char(comp)).avg(1:lengths.(char(comp)));
                total_data.(char(comp)).std = total_data.(char(comp)).std(1:lengths.(char(comp)));
%                 total_data.(char(comp)).pop = total_data.(char(comp)).data(1:lengths.(char(comp)));

                Trial_Offset_ms = size(total_data.(char(comp)).avg, 2)*(round(ms*Bins));
                x_axis = linspace(0 ,Trial_Offset_ms,size(total_data.(char(comp)).avg, 2))-PreEventNumber;

                x = x_axis';
                
                d  =total_data.(char(comp)).avg;
                sd = total_data.(char(comp)).std;      
                
                data4analyze = total_data.(char(comp)).data(1:lengths.(char(comp)));

                n  = [];
                ci = [];
                for i=1:size(data4analyze, 2)
                    n(i) = sum(~isnan(total_data.(char(comp)).data(:,i)));
                end
                
                se = sd./(n'.^0.5)';
                t = d./se;
                for i=1:size(data4analyze, 2)
                    ci(:, i) = tinv(0.975, n(i)-1)*se(i);
                end
                
%                 N = size(total_data.(char(comp)).data, 1);
%                 se = sd./(N^0.5)';
%                 t = d./se;
%                 ci = tinv(0.975, N-1)*se;
                
                hold on;
                
                ci_data = total_data.(char(comp)).data';
                ci_data = ci_data(1:length(x), :);
               
                q = 0.5;
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

                        [~,idxMin] = nanmin(ciWidth);
                        HDImin(s, :) = sortedVec(idxMin);
                        HDImax(s, :) = sortedVec(idxMin + ciIdx);
                    catch
                        HDImin(s, :) = nan;
                        HDImax(s, :) = nan;
                    end
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

%                 h = fill([x_axis';flipud(x_axis')],[d'-scattering_data';flipud(d'+scattering_data')], cmap3(comp_id,:), 'LineStyle', '--', 'EdgeColor', cmap2(comp_id,:));
%                 set(h,'facealpha',.2)

                h = fill([x_axis';flipud(x_axis')],[d'-scattering_data';flipud(d'+scattering_data')], cmap3(comp_id,:), 'LineStyle', '-', 'EdgeColor', cmap2(comp_id,:), 'Parent', fig);
                set(h,'facealpha',.5)
            end
            
            
            for i = 1:size(comp_names, 1)
                comp = char(comp_names(i));
                total_data.(char(comp)).avg = nanmean(total_data.(char(comp)).data, 1);
                total_data.(char(comp)).std = nanstd(total_data.(char(comp)).data, 1);
                % remove the last bin
                total_data.(char(comp)).avg = total_data.(char(comp)).avg(1:end);
                total_data.(char(comp)).std = total_data.(char(comp)).std(1:end);

                event_names = fieldnames(total_data.(char(comp)).events);
                for j=1:size(event_names, 1)
                    evnt = char(event_names(j));
                    total_data.(char(comp)).avg_events.(evnt) = nanmean(total_data.(char(comp)).events.(evnt));
                end
                total_data.(char(comp)).avg = total_data.(char(comp)).avg(1:lengths.(char(comp)));
                total_data.(char(comp)).std = total_data.(char(comp)).std(1:lengths.(char(comp)));

                Trial_Offset_ms = size(total_data.(char(comp)).avg, 2)*(round(ms*Bins));
                x_axis = linspace(0 ,Trial_Offset_ms,size(total_data.(char(comp)).avg, 2))-PreEventNumber;

                x = x_axis';
                
                hold on;
                plot(x_axis, total_data.(char(comp)).avg, 'LineWidth', 3, 'Color', cmap2(i,:), 'Parent', fig);
            end
%             comp1 = char(comp_names(1));
%             comp2 = char(comp_names(2));
%             total_data.diff.avg=total_data.(char(comp2)).avg(1:5500)-total_data.(char(comp1)).avg(1:5500);
%             plot(x_axis(1:5500), total_data.diff.avg, 'LineWidth', 3, 'Color', cmap2(3,:), 'Parent', fig);

            total_data.x_axis = x_axis;
            %% plot events
            for i = 1:size(comp_names, 1)
                comp = char(comp_names(i));
                Trial_Onset = 1;
                if isfield(total_data.(char(comp)).avg_events, ['event_' char(from)])
                    Trial_Onset = total_data.(char(comp)).avg_events.(['event_' char(from)]);
                end
                event_names = fieldnames(total_data.(char(comp)).events);

                for j=1:size(event_names, 1)-1
                    evnt = char(event_names(j));
                    total_data.(char(comp)).avg_events.(evnt) = mean(total_data.(char(comp)).events.(evnt));
                    x = total_data.(char(comp)).avg_events.(evnt);
                    x = x-Trial_Onset;
                    y = get(fig,'ylim');
                    plot([x x], y, 'Color', cmap2(i, :), 'LineStyle', '--', 'Parent', fig);
                end
            end

            y_unit_label   = 'arbitrary units';
            relative_label = '';
            if strcmp(relative, 'difference') || strcmp(relative, 'percentage')
                relative_label = 'Relative ';
            end
            if strcmp(Method,'mm')
                y_unit_label = 'mm';
            end
            if strcmp(Method,'Z score')
                y_unit_label = 'Z score';
            end
            if strcmp(relative, 'percentage')
                y_unit_label = '%';
            end
            xlabel_text = 'Time [ms]';
            xlabel(xlabel_text, 'FontWeight','bold');
            ylabel_text = [relative_label 'Pupil Size [' y_unit_label ']'];
            ylabel(ylabel_text, 'FontWeight','bold');
            set(gca,'FontWeight','bold');
            set(gca,'box','on');            

            title_text = 'All Participants';
            project_name = strrep(configuration.paths.project_name, '_', ' ');

            if( ~strcmp(project_name, ''))
                title_text = strcat(title_text, ' (', project_name, ')');
            end
            title(title_text);
            comp_names_fixed = cellfun(@(x) x(3:end), comp_names, 'UniformOutput', false);
            comp_names_fixed = strrep(strrep(comp_names_fixed, '_x_', ' & '),'_',' ');
            legend(fig, char(comp_names_fixed), 'Location', 'Best');
            set(fig, 'LineWidth', 2);
            xtickformat(fig, '%,.4g');

            total_data.xlabel = xlabel_text;
            total_data.ylabel = ylabel_text;
            total_data.title  = title_text;
            xlim(fig, [x_axis(1) x_axis(end)])            
        end
    end
end