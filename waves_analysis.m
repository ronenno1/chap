classdef waves_analysis
    methods(Static)        
        function run_waves_analysis()
            [data, folder, files, condition2analyze, condition_names, log_scale] = waves_analysis.load_files;
            if isempty(data)
                return;
            end
            if isempty(files)
                disp('Error: files cannot be found!');
                return;
            end

            [ps, fs] = waves_analysis.parse_all_data(data, [folder filesep, 'mat'], files, condition2analyze);

            close(findobj('type', 'figure', 'name', 'Waves analysis'))

            fig = figure('Name', 'Waves analysis');
            hold on;

            waves_analysis.print_figure(ps, fs, condition_names, log_scale);
            waves_analysis.save_data_table(ps, fs, condition_names, folder, files);
%%% Possible analysis... Note, the frequencies might be different for each condition %%%%
%             if size(ps, 1) > 1
%                 for cond_1 = 1: size(ps, 1)
%                     for cond_2 = cond_1+1: size(ps, 1)
%                         d1 = ps{cond_1, :};
%                         d2 = ps{cond_2, :};
%                         f1 = fs{cond_1, :};
%                         f2 = fs{cond_2, :};
% 
%                         num_of_frequencies = min(length(ps{1, :}'), length(ps{2, :}'));
%                         for frequency = 1:num_of_frequencies
%                             data1 = d1(:, frequency);
%                             data2 = d2(:, frequency);
%                             
%                             [t, bf, n, sd, pes, pvalue]      = stat.ttest_and_bf(data1, data2);
%                             stat_data.fs(frequency, :)       = 0.5*(f1(1, frequency) + f2(1, frequency));
%                             stat_data.se(frequency, :)       = sd./(n.^0.5);
%                             stat_data.pes(frequency, :)      = pes;
%                             stat_data.tValues(frequency, :)  = t;
%                             stat_data.pValues(frequency, :)  = pvalue;
%                             stat_data.BFs(frequency, :)      = bf;
%                         end
%                         writetable(struct2table(stat_data), [folder, 'waves', filesep, condition_names{cond_1}, ' vs. ', condition_names{cond_2}, '.csv'])
%                     end
%                 end
%             end
            file_name = [folder, 'waves', filesep, 'output'];
            if log_scale
                file_name = [file_name, '_log_scale'];
            end
            savefig(fig, file_name);
            print(fig, file_name, '-dpng', '-r300');
            close(findobj('type', 'figure', 'name', 'Waves analysis'))
        end
        
        function save_data_table(ps, fs, condition_names, folder, files)
            path2save = [folder, filesep, 'waves'];
            if ~exist([folder, filesep, 'waves'], 'file')
                mkdir(path2save);
            end

            for cond_id = 1:length(condition_names)
                cond_data = ps{cond_id};
                files     = cellfun(@(x) x(1:end-4), files, 'UniformOutput', false);

                files_data = array2table([{'frequency(Hz)'}; files], 'VariableNames', {'id'});
                freq_data  = array2table([fs{cond_id}(2, 2:end); cond_data(:, 2:end)]);
                
                writetable([files_data, freq_data], [path2save, filesep, condition_names{cond_id}, '.csv']);
            end
                
        end
        
        function [data, base_path, files, condition2analyze, condition_names, log_scale] = load_files()
            data              = [];
            files             = [];

            condition2analyze = [];
            condition_names   = [];
            log_scale         = 0;
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
            full_condition_names = data.configuration.comp_names;
            condition_names = cellfun(@(x) x(3:end), full_condition_names, 'UniformOutput', false);
            condition_names = strrep(strrep(condition_names, '_x_', ' & '), '_', ' ');
            required_condition        = listdlg('PromptString',...
                                       'Select condition:',...
                                       'ListSize',[200, 100],...
                                       'SelectionMode','multiple',...
                                       'ListString', condition_names);    
            if isempty(required_condition)
                condition2analyze = full_condition_names;
            else
                condition2analyze = full_condition_names(required_condition);
            end

            condition_names = cellfun(@(x) x(3:end), condition2analyze, 'UniformOutput', false);
            condition_names = strrep(strrep(condition_names, '_x_', ' & '), '_', ' ');
            
            folder = [base_path, filesep, 'mat'];
            files = dir([folder, filesep, '*.mat']);
            files = {files.name}';
            
            scale     = listdlg('PromptString',...
                           'Select scale:',...
                           'ListSize',[200, 100],...
                           'SelectionMode','single',...
                           'ListString', [{'Linear scale'}; {'Logarithmic scale'}]);    
            log_scale = scale==2;

        end

        function [ps, fs] = parse_all_data(data, folder, files, condition2analyze)
            ps = cell(length(condition2analyze), 1);
            fs = cell(length(condition2analyze), 1);

            for condition = 1:length(condition2analyze)
                min_duration = inf;
                for file = 1:length(files)
                    data_file    = load([folder, filesep, files{file}]);
                    cond_data    = data_file.ploted_data.(condition2analyze{condition}).cuted_data';
                    min_duration = min(min_duration, size(cond_data, 2));
                end
                
                num_of_freq = ceil(min_duration/2);

                if ~mod(min_duration, 2) % have to be odd
                    min_duration = min_duration-1;
                end

                file_ps = zeros(length(files), num_of_freq);
                file_fs = zeros(length(files), num_of_freq);
                
                for file = 1:length(files)
                    data_file = load([folder, filesep, files{file}]);
                    cond_data = data_file.ploted_data.(condition2analyze{condition}).cuted_data';
                    local_ps = zeros(size(cond_data, 1), num_of_freq);
                    local_fs = zeros(size(cond_data, 1), num_of_freq);
                    for trial = 1:size(cond_data, 1)
                        trial_data  = cond_data(trial, 1:min_duration);
                        [p, f]      = waves_analysis.parse_data(trial_data, data.configuration.rate);
                        local_ps(trial, :) = p';
                        local_fs(trial, :) = f;
                    end
                    
                    file_ps(file, :) = mean(local_ps, 'omitnan')';
                    file_fs(file, :) = mean(local_fs, 'omitnan')';
                end
                ps{condition} = file_ps;
                fs{condition} = file_fs;
            end
        end

        
        function [ps, fs] = parse_file(data, condition2analyze)
            ps = cell(length(condition2analyze));
            fs = cell(length(condition2analyze));

            for condition = 1:length(condition2analyze)
                [p, f] = waves_analysis.parse_data(data, condition2analyze{condition});
                ps{condition} = p;
                fs{condition} = f;
            end
        end

        function [p, f] = parse_data(pupil_data, Fs)
            if ~mod(length(pupil_data), 2)
                pupil_data = pupil_data(:, 1:end-1);
            end
            [p, f] = ffts.do_fft(pupil_data, Fs);
        end
        
        
        function print_figure(ps, fs, condition_names, log_scale)
            if ~exist('log_scale', 'var')
                log_scale = false;
            end

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
                [x, y, se] = waves_analysis.calculate_params(ps{condition}, fs{condition}, log_scale);
                h = fill([x fliplr(x)], [y-se  fliplr(y+se)], cmap3(condition,:),  'LineStyle', ':', 'EdgeColor', cmap2(condition,:), 'HandleVisibility', 'off');
                set(h,'facealpha',.5);
            end
    
            for condition = 1:length(condition_names)
                [x, y, ~] = waves_analysis.calculate_params(ps{condition}, fs{condition}, log_scale);
                plot(x, y, 'Color', cmap2(condition, :), 'DisplayName', condition_names{condition});
            end
            if log_scale
                xlabel('log_{10}(f) [Hz]')
                ylabel('|log_{10}(P(f))|')
                xlim([-inf, 0]);
            else
                xlabel('f [Hz]')
                ylabel('P(f)')
                xlim([-inf, 1]);
            end
            legend();      
        end
        
        function [x, y, se] = calculate_params(p, f, log_scale)
            n  = size(p, 1);
            if log_scale
                x  = mean(log10(f));
                y  = mean(log10(p));
                se = std(log10(p))./(n.^0.5)';
            else
                x  = mean(f);
                y  = mean(p);
                se = std(p)./(n.^0.5)';
            end
            x  = x(2:end);
            y  = y(2:end);
            se = se(2:end);
        end
    end
end