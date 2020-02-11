function chap()
    close all;
    clc
    toolboxes = ver;
    toolbox_names = {toolboxes.Name};
    if(isempty(find(~cellfun(@isempty, strfind(toolbox_names,'Curve Fitting Toolbox')),1)))
        disp('It seems that you do not have the Curve Fitting Toolbox. It is highly recommended to install it.');
    end
    Figure_h = gui_lib.create_figure('CHAP - Main', 'images/chap_1.jpg', [0 0 550 440]);%create the initial figure
    add_main_buttons(Figure_h)  %create for buttons and present them in the figure
end

function add_main_buttons(Figure_h)
    log = text(10, 420, '', 'HorizontalAlignment', 'left', 'Color', [1 1 1], 'FontWeight', 'bold');

    gui_lib.uicontrol_button(Figure_h, [22 332 195 40], 'New project', {@read_raw_file, log});
    gui_lib.uicontrol_button(Figure_h, [22 280 195 40], 'Open existing project', {@read_chap, log});
    gui_lib.uicontrol_button(Figure_h, [22 228 195 40], 'Between groups analysis', {@idttest, log});

%     gui_lib.uicontrol_button(Figure_h, [228 332 195 40], 'Convert EDF to Matlab', {@quickEdf2mat, log});    
    gui_lib.uicontrol_button(Figure_h, [228 332 195 40], 'About', @about);    

    gui_lib.uicontrol_button(Figure_h, [22 65 125 30], 'Exit', @exit);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   primary functions                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_raw_file(src, evnt, log) % read data file and convert it into mat file
    data = read_data.read_file(log);
    if isempty(data)
        return;
    end
    guidata(src, data);
    process_window(src, evnt, log);
end

function about(~, ~)
    url = 'http://in.bgu.ac.il/en/Labs/CNL/chap';
    try
        web(url, '-browser');
    catch
        web(url);
    end
end


% function quickEdf2mat(src, evnt, log) %read EDF file and convert it to mat file
%     start                   = tic;
%     [edf_file_name, edf_path_name]  = uigetfile('*.edf','Select the EDF file');
%     full_edf_name                   = [edf_path_name edf_file_name];
%     if ~edf_path_name
%         return;
%     end
% 
%     [mat_file_name ,mat_path_name] = uiputfile('*.mat', 'Select file to save');
%     full_mat_name               = [mat_path_name mat_file_name];
%     if(~mat_file_name)
%         return;
%     end
%     
%     print_log(['"' strrep(edf_file_name,  '_', '\_') '" start convert: ' full_edf_name], log);    
%     edf2simpleMat(full_edf_name, full_mat_name, log);
%     print_log(['"' strrep(edf_file_name,  '_', '\_') '" successfully convert! ' num2str(toc(start)) ' seconds'], log);    
% end

function idttest(src, evnt, log) % read edf file and convert it to mat file
    idttest.run_idttest;
end

function read_chap(src, evnt, log) % read edf file and convert it to mat file
    start = tic;
    [chap_file_name, chap_path_name]   = uigetfile('*.chp','Select the Chap file');
    full_chap_name                   = [chap_path_name chap_file_name];
    if(~chap_path_name)
        return;
    end
    cd(chap_path_name)
    
    file      = load(full_chap_name,'-mat');
    chap_data = file.data;
    data      = read_data.parse_data(chap_data, log);
    
    guidata(src, data);
    print_log([strrep(chap_file_name, '_', '\_') ' successfully load! ' num2str(toc(start)) ' seconds'], log);    
    process_window(src, evnt, log);
end

function exit(~, ~) % exit
    close;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    utility functions                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function process_window (src, ~, log)
    data = guidata(src);
    if isempty(data)
        print_log('data not found', log);    
        return
    end
    conds_data = data.conds_data;
    data_mat   = data.data_mat;
    vars       = data_mat.Properties.VariableNames();
    
    %% init figure
    figure = gui_lib.create_figure('CHAP - Data Pre-processing', 'images/chap_2.jpg', [0 0 900 750]); %create gui figure
    
    data = view_configuration(data);

    handler = view_user_vars(figure, conds_data, vars);
    
    if(isempty(handler))
        print_log('variables not found', log);    
        close(findobj('type','figure','name','CHAP - Data Pre-Processing'))
        return;
    end
    
    log = text(10, 730, '', 'HorizontalAlignment', 'left', 'Color', [1 1 1], 'FontWeight', 'bold');

    gui_lib.uicontrol_button(figure, [32 40 125 30], 'Exit', @exit);
    gui_lib.uicontrol_button(figure, [162 40 125 30],'Continue', {@start_analyze, data, handler, log}); 
end

function start_analyze(src, ~, data, handler, log)
    conds_data         = data.conds_data;
    data.var_handler   = save_vars(handler, conds_data);
    
    data.configuration = Get_processing_vals(data.configuration);
    
    data = analyze_data(data, data.configuration.min_trials_val, data.configuration.ZOutliers_val, data.configuration.ZeroshNumber_val, [], log);
    if isempty(data) || ~isfield(data, 'cond_mat_cleaned')
        print_log('Error: please select variables', log);    
        return
    end
    
    cond_mat_cleaned   = data.cond_mat_cleaned;
    cond_events        = data.cond_events;
    data               = show_analyze_window(src, cond_mat_cleaned, cond_events, data);
    guidata(src, data);
end

function do_analyze(data, fig, src, log, log_a)

    data.configuration = get_analyze_vals(data.configuration);
    analyzed_data      = reload_graph(data, data.configuration, fig);
    error_msg = 'Wrong range (too many samples before the first event)';
    if (isfield(analyzed_data, 'wrong_range'))
        print_log('Error: ', log_a);   
        print_log(error_msg, log);    
    elseif strfind(log.String, error_msg)
        print_log('', log_a);   
        print_log('', log);    
    end
    data.analyzed_data = analyzed_data;
    guidata(src, data);
end


function data = analyze_data(data, min_trials, z_outliers, zeros_outliers, comp_names, log)
    [data_mat, blinks_data] = fix_all_blinks(data, z_outliers, zeros_outliers, log);
    data_mat.pupil(data_mat.pupil==0) = nan;    
    data.blinks_data = blinks_data;
    data.data_mean  = nanmean(data_mat.pupil);
    data.data_std   = nanstd(data_mat.pupil);
    data.min_trials = min_trials;

    var_handler  = data.var_handler; % data about handlers - what user choose

    conds       = var_handler.conditions;
    event_names = var_handler.events;
    data.configuration.event_names = [event_names 'Trial_Offset'];
    if(isempty(conds))
        return;
    end
    cond_names  = fieldnames(conds);
    comp_nums   = analyzer.make_comps(conds); % return index variables for all the posible comparations
    print_log('Start parsing levels', log);  
    data = analyzer.analyze_all(data, data_mat, cond_names, comp_nums, conds, event_names, comp_names, log);
end

function [data_mat, blinks_data] = fix_all_blinks(data, z_outliers, zeros_outliers, log)
    data_mat = data.data_mat;
    view_trials        = gui_lib.get_checkbox_val(data.configuration.view_trials);
    interpolation_type = gui_lib.get_popupmenu_val(data.configuration.interpolation_type);
    
    linear_interpolation = false;
    if strcmp(interpolation_type, 'Linear Interpolation')
        linear_interpolation = true;
    end
    if(data_mat.pupil_size(1, 1)<0) % eyelink support only one eye
        eye = 2; %right
    else
        eye = 1;
    end
    data_mat.pupil = zeros(size(data_mat,1),1);
    
    
    gradient_crit = 4;
    nan_data = data_mat.pupil_size(:,eye);
    nan_data(nan_data==0)=nan; 
%     total_mean = nanmean(nan_data);
%     total_std = nanstd(nan_data);
    diff_data = diff(nan_data);
    diff_mean = nanmean(diff_data);
    diff_std  = nanstd(diff_data);
    
    gradient    = diff_mean+gradient_crit*diff_std;
        
    data_mat(data_mat.trial_id<0, :) = [];
    
%     all_pupil_data = data_mat.pupil_size(:,eye); 
%     [data_mat.pupil_size(:,eye), ~]    = fix_blinks2(all_pupil_data, z_outliers, zeros_outliers, data.rate, linear_interpolation, gradient, true);

    [trials, trial_ids] = unique(data_mat.trial_id);
    
    percentages = zeros(1, 10);
    print_log('Start blinks fixing', log);
    blinks_data = cell(size(trials));

    for trial=1:size(trials, 1)

        % logger
        percentage = round(100*(trial/size(trials,1)));
        if percentage>0 && ~mod(percentage, 10) && percentages(percentage/10)==0 
            print_log(['Blinks fixed: ' num2str(percentage) '%'], log);    
            percentages(percentage/10) = 1;
        end
        % logger end
        Trial_Offset = size(data_mat, 1);
        if trial < size(trials, 1)
            Trial_Offset = trial_ids(trial+1)-1;
        end
        Trial_Onset = trial_ids(trial);
        trial_data = data_mat(Trial_Onset:Trial_Offset, :);
        trial2show = 0;
        if view_trials
            trial2show = trial;
        end
        if strcmp(interpolation_type, 'Without Interpolation')
            data_mat.pupil(Trial_Onset:Trial_Offset) = trial_data.pupil_size(:,eye);
            continue;
        end
        pupil_data = trial_data.pupil_size(:,eye);
        [pupil_data, blinks_data_positions]       = fix_blinks2(trial_data.pupil_size(:,eye), z_outliers, zeros_outliers, data.rate, linear_interpolation, gradient, trial2show);
        data_mat.pupil(Trial_Onset:Trial_Offset)  = pupil_data;
        blinks_data{trial} = blinks_data_positions;
    end
end

function change_comps(src, cond_mat_cleaned, cond_events)
    data = guidata(src);
    close(findobj('type','figure','name','CHAP - Advanced Processing Options and Time-Course Visualization'))
    close(findobj('type','figure','name','CHAP - Statistical Analysis'))

    show_analyze_window(src, cond_mat_cleaned, cond_events, data);
end

function [ploted_data, events] = reload_graph(data, configuration, fig)
    [ploted_data, events] = ploter2.draw_graph(data.cond_mat_data, data.cond_events_data, data.cond_blinks, configuration.comp_names, data.rate, fig, data.data_mean, data.data_std, configuration);
    title(strrep(data.file_name, '_', '\_'));        
end

function [figure, log, log_a, fig, analyze_bot] = analyze_window()
    figure = gui_lib.create_figure('CHAP - Advanced Processing Options and Time-Course Visualization', 'images/chap_3.jpg', [0 0 900 750]);
    gui_lib.uicontrol_title('CONFIGURATION', 130, 70);
    gui_lib.uicontrol_title('TIME COURSE', 564, 70);
    gui_lib.uicontrol_title('DESCRIPTIVE INFORMATION', 358, 493);
    log_a = text(70, 730, '', 'HorizontalAlignment', 'right', 'Color', [1 1 1], 'FontWeight', 'bold');
    log   = text(70, 730, '', 'HorizontalAlignment', 'left', 'Color', [1 1 1], 'FontWeight', 'bold');
    hp	  = uipanel('Units','Pixels', 'BackgroundColor',[1 1 1], 'Position',[362 325 500 340]);    
    fig   = subplot(1, 1, 1, 'Parent', hp);
    analyze_bot = gui_lib.uicontrol_button(figure, [32 281 300 30], 'Analyze');
end

function [figure, log, log_a, fig, analyze_bot] = statistical_window()
    close(findobj('type','figure','name','CHAP - Statistical Analysis'))

    figure = gui_lib.create_figure('CHAP - Statistical Analysis', 'images/chap_3.jpg', [0 0 900 750]);
    gui_lib.uicontrol_title('CONFIGURATION', 130, 70);
    gui_lib.uicontrol_title('RESULTS', 578, 70);
    gui_lib.uicontrol_title('RESULTS', 418, 493);
    log_a = text(70, 730, '', 'HorizontalAlignment', 'right', 'Color', [1 1 1], 'FontWeight', 'bold');
    log   = text(70, 730, '', 'HorizontalAlignment', 'left', 'Color', [1 1 1], 'FontWeight', 'bold');
    hp	  = uipanel('Units','Pixels', 'BackgroundColor',[1 1 1], 'Position',[362 325 500 340]);    
    fig   = subplot(1, 1, 1, 'Parent', hp);
    analyze_bot = gui_lib.uicontrol_button(figure, [32 281 300 30], 'Analyze');

    gui_lib.uicontrol_button(figure, [32 40 125 30], 'Exit', @exit);
end


function show_statistical_window(~, statistical_data, total_data, comp_names)
         
    [figure, ~, ~, fig, analyze_bot]  =  statistical_window();
       
    statistical_data = show_statistical_vars(statistical_data);
    
    tooltip = '<html><i>Save the data of the participant (mat file)';
    next = gui_lib.uicontrol_button(figure, [714 281 155 30], 'Next', tooltip);   
    back = gui_lib.uicontrol_button(figure, [534 281 155 30], 'Back', tooltip);  
    
    tooltip = '<html><i>Save the figure as image (png file) or as matlab figure (fig file)';
    
    gui_lib.uicontrol_button(figure, [354 281 155 30], 'Save figure', @(~,~)output.save_figure2(fig), tooltip);
    stat.plot_statistics(0,0, fig, figure, statistical_data, total_data, 1, 0, next, back);
    set(analyze_bot, 'callback', {@stat.calc_and_plot, fig, figure, total_data, statistical_data, comp_names, next, back, statistical_data.configuration});
end


function data = show_analyze_window(src, cond_mat_cleaned, cond_events, data)
    comp_names       = fieldnames(cond_mat_cleaned);

    comp_names_fixed = cellfun(@(x) x(3:end), comp_names, 'UniformOutput', false);
    header           = strrep(strrep(comp_names_fixed, '_x_', ' & '),'_',' ');
    DataValue        = listdlg('PromptString',...
                               'Select conditions:',...
                               'ListSize',[300 300],...
                               'SelectionMode','multiple',...
                               'ListString', header);    
    if isempty(DataValue)
        return;
    end
    
    comp_names = comp_names(DataValue);
    [figure, log, log_a, fig, analyze_bot]  =  analyze_window();
    
    [ploted_data, ~]    = ploter2.draw_graph(data.cond_mat_data, data.cond_events_data, data.cond_blinks, comp_names, data.rate, fig, data.data_mean, data.data_std, data.configuration);

    data.analyzed_data = ploted_data;
    
    data = show_analyze_vars(data);
    
    data.configuration.comp_names = comp_names;

    outliers                      = data.outliers;
    valid_trials                  = data.valid_trials;
   
    set(analyze_bot, 'callback', @(~,~)do_analyze(data, fig, src, log, log_a)); 

    num_of_events = size(fieldnames(ploted_data.(char(comp_names(1))).events),1);
    data_table    = zeros(size(comp_names,1), num_of_events + 2);
    
    data.fail = false;
    for i=1:size(comp_names, 1)
        num_of_valid_trials = valid_trials.(char(comp_names(i)));
        data_table(i,:) = [outliers.(char(comp_names(i))) num_of_valid_trials int32(struct2array(ploted_data.(char(comp_names(i))).events))];
        if(num_of_valid_trials<data.min_trials)
            print_log('Error: ', log_a);   
            print_log('not enough trials', log);    
            data.fail = true;
        end
    end

    view_analyze_results(ploted_data, comp_names, data_table, src, cond_mat_cleaned, cond_events, figure, fig, log, log_a);
    do_analyze(data, fig, src, log, log_a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      group functions                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function do_group(src, fig, log, log_a, res_table)

    print_log('', log_a);    
    print_log('', log);    

    data = guidata(src);
    do_analyze(data, fig, src, log, log_a);
    data = guidata(src);
    [paths, files, overwrite] = grouper.init_folders(data.file_type);
    if(isempty(paths))
        return;
    end
    if(isempty(files))
        print_log('Error:', log_a);    
        print_log('files not found', log);    
        return;
    end
    configuration = data.configuration;
    comp_names    = configuration.comp_names;
    configuration.paths = paths;
    process_files(paths, files, data, fig, log, log_a, res_table, overwrite);
   
    if size(files, 1)<2
        print_log('Error:', log_a);    
        print_log('You need at least 2 files for group analysis', log);    
        return;
    end

    total_data = grouper.do_avaraging(paths, comp_names);
    if (isempty(total_data) || size(total_data.(char(comp_names(1))).data, 1)<2)
        
        print_log('Error:', log_a);    
        print_log('There are no enough valid files', log);    
        return;
    end

    print_log('Combine all the data...', log);    
    [total_data, comp_names_fixed] = grouper.do_plot(total_data, configuration, data.rate, fig);
    
    header  = fieldnames(total_data.(char(comp_names(1))).events)';
    header = cellfun(@(x) iff(strcmp(x(1:6),'event_'),  x(7:end) ,x), header, 'UniformOutput', false);
    header = strrep(header, '_', ' ');
 
    num_of_events = size(header, 2);
    data_table    = zeros(size(comp_names,1), num_of_events);
    
    for i=1:size(comp_names, 1)
        data_table(i,:) = int32(struct2array(total_data.(char(comp_names{i})).avg_events));
    end
    set(res_table, 'Data', data_table);
    set(res_table, 'ColumnName', header');
    total_data.paths = paths;
    configuration.rate = data.rate;
    total_data.configuration = configuration;
    
    output.save_figure(fig, comp_names_fixed, [paths.png_output_folder_name filesep 'total.png']);
    output.save_figure(fig, comp_names_fixed, [paths.fig_output_folder_name filesep 'total.fig']);
  
    save([paths.project_output_folder_name filesep 'total'], 'total_data');
    properties.Z_outliers = total_data.configuration.ZOutliers_val;
    properties.missing_values = total_data.configuration.ZeroshNumber_val;
    properties.min_trials = total_data.configuration.min_trials_val;
    properties.interpolation_type = gui_lib.get_popupmenu_val(total_data.configuration.interpolation_type);
    properties.bins = total_data.configuration.BinsNumber_val;
    properties.from = total_data.configuration.from_val;
    properties.to = total_data.configuration.to_val;
    properties.pre_event_ms = total_data.configuration.PreEventNumber_val;
    properties.relative = total_data.configuration.relative_val;
    properties.method = total_data.configuration.Method_val;
    properties.baseline = total_data.configuration.baseline_val;
    properties.scattering = total_data.configuration.scattering_val;
   
    event_names  = cellfun(@(e) {e(7:end)}, total_data.configuration.event_names, 'UniformOutput', false);
    properties.event = event_names(1:end-1);
    contrast_names_fixed   = strrep(strrep(strrep(total_data.configuration.comp_names_val', 'c_', ''), '_x_', ' & '),'_',' ');

    properties.var = contrast_names_fixed;

    properties_table = struct2table(properties);
    writetable(properties_table, [paths.project_output_folder_name filesep 'properties.csv']);  

    if (size(files, 1)>1 && length(comp_names)>1)
        statistic_data = stat.do_stat(total_data, comp_names);
        show_statistical_window(src, statistic_data, total_data, comp_names);
    end
    print_log('', log_a);    
    print_log('Done!', log);    
end

function process_files(paths, files, data, fig, log, log_a, res_table, overwrite)
    
    %% load data from previous analysis (in case it will be required to use it) 
    if ~overwrite
        if exist(strcat(paths.behave_output_folder_name, filesep, 'time-course_data.csv'), 'file')
            old_behave_table = readtable(strcat(paths.behave_output_folder_name, filesep, 'time-course_data.csv'));
        end
        if exist(strcat(paths.behave_output_folder_name, filesep, 'outliers.csv'), 'file')
            old_fail_behave_table = readtable(strcat(paths.behave_output_folder_name, filesep, 'outliers.csv'));
        end
        
        if exist(strcat(paths.csv_output_folder_name, filesep, 'trials_data.csv'), 'file')
            old_var_data_table = readtable(strcat(paths.csv_output_folder_name, filesep, 'trials_data.csv'));
        end
        
        if exist(strcat(paths.csv_output_folder_name, filesep, 'outliers.csv'), 'file')
            old_outliers = readtable(strcat(paths.csv_output_folder_name, filesep, 'outliers.csv'));
        end
        
        if exist(strcat(paths.csv_output_folder_name, filesep, 'valid_trials.csv'), 'file')
            old_valid_trials = readtable(strcat(paths.csv_output_folder_name, filesep, 'valid_trials.csv'));
        end
    end
    
    configuration = data.configuration;

    num_of_files  = size(files, 1);
    if overwrite || ~exist(strcat(paths.csv_output_folder_name, filesep, 'time-course_data.csv'), 'file')
        output.save_csv_header(configuration.comp_names, [paths.csv_output_folder_name filesep 'time-course_data.csv']);
    end
    comp_names        = configuration.comp_names;
    var_data_table    = [];
    behave_table      = [];
    fail_behave_table = [];
    updated = false;
    for i = 1:num_of_files
        print_log([num2str(i) '/' num2str(num_of_files) ':'], log_a);
        print_log(['Loading: ' char(files(i)) ], log);
        [~, file_name, ~] = fileparts(char(files(i)));
        
        skip_me = ~overwrite && (exist([paths.mat_output_folder_name filesep file_name '.mat'], 'file') || exist([paths.mat_output_folder_name_err filesep file_name '.mat'], 'file'));
    
        if skip_me && exist('old_outliers', 'var')
            participant_id = strcmp(old_outliers.participant_id, file_name);
            participants.outliers(i) = table2struct(old_outliers(participant_id, :));
        end
        
        if skip_me && exist('old_valid_trials', 'var')
            participant_id = strcmp(old_valid_trials.participant_id, file_name);
            participants.valid_trials(i) = table2struct(old_valid_trials(participant_id, :));
        end

        if skip_me
            continue;
        end
        updated = true;
        %% get data from file        
        single_data             = process_file([paths.output_folder_name filesep char(files(i))], paths.chap_output_folder_name, log, data.events2, data.vars2);
        single_data.var_handler = data.var_handler;
        single_data.configuration = configuration;
        %% ignore
        if(~isfield(single_data, 'data_mat'))
            continue;
        end
        single_data             = analyze_data(single_data, configuration.min_trials_val, configuration.ZOutliers_val, configuration.ZeroshNumber_val, comp_names, log);
        if isempty(single_data)
            print_log(['Error: wrong file: ' char(files(i)) ': '], log);    
            return;
        end
        if(single_data.rate~=data.rate)
            print_log(['Error: wrong smapling rate: ' char(files(i)) ': '], log);    
            continue;
        end
        skip = false;
        if(~isfield(single_data, 'cond_mat_data'))
           print_log(['Error: no data: ' char(files(i)) ': '], log);    
           continue;
        end

        for comp = 1:size(comp_names, 1)
            if(~isfield(single_data.cond_mat_data, char(comp_names(comp))))
               print_log(['Error: no data for ' char(comp_names(comp))  ': ' char(files(i)) ': '], log);    
                skip = true;
             end
        end
        if skip
            continue;
        end
        cond_ids            = single_data.cond_ids;     
        [ploted_data, ~]    = ploter2.draw_graph(single_data.cond_mat_data, single_data.cond_events_data, single_data.cond_blinks, comp_names, single_data.rate, fig, single_data.data_mean, single_data.data_std, data.configuration, char(files(i)));
       	if(isfield(ploted_data, 'wrong_range'))
           print_log(['Error: wrong_range: ' char(files(i)) ': '], log);    
            continue;
        end
        drawnow

        %% save the graphs
        comp_names_fixed = cellfun(@(x) x(3:end), comp_names, 'UniformOutput', false);
        comp_names_fixed = strrep(strrep(comp_names_fixed, '_x_', ' & '),'_',' ');
        printed_data = ploted_data;
        %% save the data
        comp_names = fieldnames(ploted_data);

        ploted_data.rate = data.rate;
        ploted_data.fail = false;
        for j = 1:size(comp_names,1)
            comp = char(comp_names(j));
            ploted_data.(char(comp)).outliers    = single_data.outliers.(char(comp));
            num_of_valid_trials = single_data.valid_trials.(char(comp));
            ploted_data.(char(comp)).valid_trials = single_data.valid_trials.(char(comp));
            if(num_of_valid_trials<data.min_trials)
                ploted_data.fail = true;
            end
        end
        ploted_data.blinks_data = single_data.blinks_data;
        
        valid_trials = single_data.valid_trials;
        valid_trials.participant_id = file_name;
        participants.valid_trials(i) = valid_trials;
        outliers = single_data.outliers;
        outliers.participant_id = file_name;
        participants.outliers(i) = outliers;

        %% print the data
        header = [cellstr('outliers') cellstr('Valid_Trials') fieldnames(ploted_data.(char(comp_names(1))).events)'];
        header = cellfun(@(x) iff(strcmp(x(1:6),'event_'),  x(7:end) ,x), header, 'UniformOutput', false);
        header_f = strrep(header, '_', ' ');    
        num_of_events = size(fieldnames(ploted_data.(char(comp_names(1))).events),1);
        data_table = zeros(size(comp_names,1), num_of_events + 2);
        for j=1:size(comp_names, 1)
            data_table(j,:) = [ploted_data.(char(comp_names(j))).outliers ploted_data.(char(comp_names(j))).valid_trials int32(struct2array(ploted_data.(char(comp_names(j))).events))];
        end
        set(res_table, 'Data', data_table);
        set(res_table, 'ColumnName', header_f');
        single_behave_table = [table(repmat({file_name}, size(comp_names_fixed)), 'VariableNames', {'participant_id'}),...
                        table(comp_names_fixed, 'VariableNames', {'condition'}),...
                        array2table(data_table, 'VariableNames',header)];
        if(ploted_data.fail)
            fail_behave_table = [fail_behave_table; single_behave_table];
        else
            behave_table = [behave_table; single_behave_table];
        end
        
        if(ploted_data.fail)
            output.save_figure(fig, comp_names_fixed, [paths.png_output_folder_name_err filesep file_name '.png']);
            output.save_figure(fig, comp_names_fixed, [paths.fig_output_folder_name_err filesep file_name '.fig']);
            print_log('Error:', log_a);
            print_log('not enough trials', log);    
            save([paths.mat_output_folder_name_err filesep file_name], 'ploted_data');
        else
            output.save_figure(fig, comp_names_fixed, [paths.png_output_folder_name filesep file_name '.png']);
            output.save_figure(fig, comp_names_fixed, [paths.fig_output_folder_name filesep file_name '.fig']);
            single_data.printed_data = printed_data;
            save([paths.mat_output_folder_name filesep file_name], 'ploted_data');
            output.save_csv_append(printed_data, [paths.csv_output_folder_name filesep 'time-course_data.csv'], file_name);
            single_var_data_table = output.save_csv_append4(single_data, cond_ids, ploted_data.rate);
            var_data_table = [var_data_table; single_var_data_table];
        end
    end 
    
    if ~updated 
        return;
    end
    if ~isempty(behave_table)
        if exist('old_behave_table', 'var')
            behave_table = [old_behave_table; behave_table];
        end            
        writetable(behave_table, strcat(paths.behave_output_folder_name, filesep, 'time-course_data.csv'));
    end
    
    if(~isempty(fail_behave_table))
        if exist('old_fail_behave_table', 'var')
            fail_behave_table = [old_fail_behave_table; fail_behave_table];
        end            
        writetable(fail_behave_table, strcat(paths.behave_output_folder_name, filesep, 'outliers.csv'));
    end
    
    if exist('var_data_table', 'var') && size(var_data_table, 1) > 0
        writetable(var_data_table, [paths.csv_output_folder_name filesep 'trials_data.csv']);
    end
    
    if exist('old_var_data_table', 'var')
        new_var_data_table = readtable(strcat(paths.csv_output_folder_name, filesep, 'trials_data.csv'));
        var_data_table = [old_var_data_table; new_var_data_table];
        writetable(var_data_table, [paths.csv_output_folder_name filesep 'trials_data.csv']);
    end
    
    if exist('participants', 'var')
        writetable(struct2table(participants.outliers), strcat(paths.csv_output_folder_name, filesep, 'outliers.csv'));
        writetable(struct2table(participants.valid_trials), strcat(paths.csv_output_folder_name, filesep, 'valid_trials.csv'));
    end
end

function data = process_file(full_name, output_folder_name, log, events2, vars2) % read edf file and convert it to mat file
    if (~exist('events2', 'var'))
        events2 = [];
    end
    if (~exist('vars2', 'var'))
        vars2 = [];
    end
    tic;
    
    data = [];
    [path, file_name, ext]   = fileparts(full_name);
    
    chap_file_name = [output_folder_name filesep file_name '.chp'];
    if(~exist(chap_file_name, 'file'))
        if(exist(strcat(path, filesep, file_name, '.chp'), 'file'))
            copyfile(strcat(path, filesep, file_name, '.chp'), chap_file_name);
            file = load(chap_file_name, '-mat');
            full_data_mat = file.data;    

        elseif strcmp(ext, '.edf')
            full_data_mat = edf2matlab2(full_name, output_folder_name, log, events2, vars2);
        elseif strcmp(ext, '.txt')
            full_data_mat  = etTxt2matlab(full_name, output_folder_name, log, events2, vars2);
        elseif strcmp(ext, '.asl')
            full_data_mat  = asl2matlab(full_name, output_folder_name, log, events2, vars2);
        elseif strcmp(ext, '.tbi')
            full_data_mat  = tobii2matlab(full_name, output_folder_name, log, events2, vars2);
        elseif strcmp(ext, '.plsd')
            full_data_mat  = plsd2matlab(full_name, output_folder_name, log, events2, vars2);
        elseif strcmp(ext, '.dat') 
            full_data_mat = dat2matlab(full_name, output_folder_name, log, events2, vars2);
        end
        if(isempty(full_data_mat))
            return;
        end
    else
        file = load(chap_file_name, '-mat');
        full_data_mat = file.data;    
    end
    full_data_mat       = load_chap_data(full_data_mat, log);
    data.data_mat       = full_data_mat.data;
    data.rate           = full_data_mat.rate;
    full_data_mat.var_data_table.event_Trial_Onset = zeros(size(full_data_mat.var_data_table, 1), 1);

    data.var_data_table = full_data_mat.var_data_table;
    data.file_name      = file_name;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      gui functions                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function data = view_configuration(data)
    hp = uipanel('Units', 'Pixels', 'BackgroundColor', [0.1 0.32 0.46], 'BorderWidth', 0, 'Position', [35 460 210 180]);
    title_pos_x = 10;
    edit_pos_x  = 150;
    title_pos_y = 140;
    edit_pos_y  = 150;

    gui_lib.uicontrol_title('CONFIGURATION', 75, 95);


    tooltip = '<html><b>Z outliers</b><br><i>Select the Z score value of outlier samples</html>';

    gui_lib.uicontrol_text(hp, 'Z outliers:', [title_pos_x title_pos_y 150 30]);
    ZNumber = gui_lib.uicontrol_edit(hp, [edit_pos_x edit_pos_y 50 25], tooltip, @valid_positive_number);
    set(ZNumber, 'String', configer.z_number);

    tooltip = '<html><b>Missing values</b><br><i>Select the percentage of missing values to exclude trial</html>';
    gui_lib.uicontrol_text(hp, 'Missing values(%):', [title_pos_x title_pos_y-35 150 30]);
    ZeroshNumber = gui_lib.uicontrol_edit(hp, [edit_pos_x edit_pos_y-35 50 25], tooltip, @valid_positive_number);
    set(ZeroshNumber, 'String', configer.zeros_number);
        
    
    tooltip = '<html><b>Min. Trials</b><br><i>Select the minimum number of trials for analyzing</html>';

    gui_lib.uicontrol_text(hp, 'Min. trials:', [title_pos_x title_pos_y-70 150 30]);
    min_trials = gui_lib.uicontrol_edit(hp, [edit_pos_x edit_pos_y-70 50 25], tooltip, @valid_2plus_number);
    set(min_trials, 'String', configer.min_trials);
    
    tooltip = '<html><b>Blinks currection (based on <a href = "https://doi.org/10.3758/s13428-018-01190-1">Hershman et.al., 2018</a>)</b><br> <br><i>Select how to reconstruct the missing values: <li> Linear interpolation <li>  Cubic interpolation <li> No correction</html>';

    val_options = {'Linear Interpolation', 'Cubic Interpolation', 'Without Interpolation'};

    data.configuration.interpolation_type = uicontrol(hp,'Style','popupmenu',...
        'String', val_options, 'Position',[title_pos_x title_pos_y-100 190 30], 'TooltipString', tooltip);%samples pre-event
    
    tooltip = '<html><b>View trials</b><br><i>View all trials including blinks currection - using for debug</html>';
    gui_lib.uicontrol_text(hp, 'View trials:', [title_pos_x title_pos_y-135 150 30]);
        
    data.configuration.view_trials  = gui_lib.uicontrol_smart_checkbox(hp, [edit_pos_x edit_pos_y-135 50 25], tooltip);

    data.configuration.min_trials   = min_trials;
    data.configuration.ZNumber      = ZNumber;
    data.configuration.ZeroshNumber = ZeroshNumber;    
end


function statistical_data = show_statistical_vars(statistical_data)
    hp = uipanel('Units', 'Pixels', 'BackgroundColor', [0.1 0.32 0.46], 'BorderWidth', 0, 'Position', [33 320 300 350]);
    title_pos_x = 10;
    edit_pos_x  = 195;
    edit_pos_y  = 300;

    tooltip = '<html><b>Statistical approach</b><br><i>Select the statistical approach for your analyzing</html>';
    
    gui_lib.uicontrol_text(hp, 'Statistical approach:', [title_pos_x edit_pos_y 150 30]);
    approach = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y 120 30], 'TooltipString', tooltip, 'String', [{'Bayesian'}; {'Classical'}], 'Value', 1);

    gui_lib.uicontrol_text(hp, 'To:', [title_pos_x edit_pos_y-80 50 30]);
    tooltip = '<html><b>To</b><br><i>Select last event for analyze</html>';

    range = statistical_data.configuration.range;

    to = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y-80 120 30], 'TooltipString', tooltip, 'String', range(2:end), 'Value', size(range, 2)-1);%samples pre-event


    gui_lib.uicontrol_text(hp, 'From:', [title_pos_x edit_pos_y-40 50 30]);
    tooltip = '<html><b>From</b><br><i>Select first event for analyze</html>';
    from = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y-40 120 30], 'TooltipString', tooltip, 'String', range(1:end-1), 'Value', 1, 'Callback',{@change_to, to, range});%samples pre-event
    
    
    
    gui_lib.uicontrol_text(hp, 'Parameter:', [title_pos_x edit_pos_y-120 90 30]);
    tooltip = '<html><b>Parameter</b><br><i>Select the desired parameter for analysis</html>';
    measure = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y-120 120 30], 'TooltipString', tooltip, 'String', [{'mean'}; {'peak'};{'peak latency'}; {'dip'};{'dip latency'};], 'Value', 1);
    
    
    
    gui_lib.uicontrol_text(hp, 'Output type:', [title_pos_x edit_pos_y-160 120 30]);
    tooltip = '<html><b>From</b><br><i>Select first event for analyze</html>';
    output_type = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y-160 120 30], 'TooltipString', tooltip, 'String', [{'descriptive'}; {'inference'}], 'Value', 1);

    
    gui_lib.uicontrol_text(hp, 'Presentation type:', [title_pos_x edit_pos_y-210 120 40]);
    tooltip = '<html><b>Presentation type</b><br><i>Select how to present the data (only for the Bayesian approach)</html>';
    presentation_type = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y-200 120 30], 'TooltipString', tooltip, 'String', [{'full'}; {'compact'};  {'combined'}], 'Value', 1);

    statistical_data.configuration.var.from                 = from;
    statistical_data.configuration.var.to                   = to;
    statistical_data.configuration.var.approach             = approach;
    statistical_data.configuration.var.measure              = measure;
    statistical_data.configuration.var.output_type          = output_type;
    statistical_data.configuration.var.presentation_type    = presentation_type;

end

function data = show_analyze_vars(data)
    hp = uipanel('Units', 'Pixels', 'BackgroundColor', [0.1 0.32 0.46], 'BorderWidth', 0, 'Position', [33 320 300 350]);
    title_pos_x = 10;
    edit_pos_x  = 195;
    title_pos_y = 290;
    edit_pos_y  = 300;

    tooltip = '<html><b>Bins</b><br><i>Select number of ms for each bin</html>';
    
    gui_lib.uicontrol_text(hp, 'Bins [ms]:', [title_pos_x title_pos_y 150 30]);
    BinsNumber = gui_lib.uicontrol_edit(hp, [edit_pos_x edit_pos_y 80 30], tooltip, @valid_positive_number);
    
    %% cut to events

    event_names  = data.configuration.event_names;
    if ~iscell(event_names)
        event_names = cellstr(event_names);
    end
    event_names  = ['Trial_Onset', cellfun(@(x) iff(strcmp(x(1:6),'event_'),  x(7:end) ,x), event_names, 'UniformOutput', false)]';    

    gui_lib.uicontrol_text(hp, 'To:', [title_pos_x edit_pos_y-80 50 30]);
    tooltip = '<html><b>To</b><br><i>Select last event for analyze</html>';

    to = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y-80 120 30], 'TooltipString', tooltip);%samples pre-event

    gui_lib.uicontrol_text(hp, 'From:', [title_pos_x edit_pos_y-40 50 30]);
    tooltip = '<html><b>From</b><br><i>Select first event for analyze</html>';
    from = uicontrol(hp,'Style','popupmenu',...
        'Position',[edit_pos_x-40 edit_pos_y-40 120 30],'Callback',{@change_to, to, event_names}, 'TooltipString', tooltip);%samples pre-event
    
    set(from, 'String', event_names(1:size(event_names, 1)-1));
    set(to, 'String', event_names(2:end));
    set(to, 'Value', size(event_names, 1)-1);
 
    gui_lib.uicontrol_text(hp, 'Ms. before event:', [title_pos_x title_pos_y-120 150 30]);
    tooltip = '<html><b>Samples before event</b><br><i>Select number of ms before events';
    PreEventNumber = gui_lib.uicontrol_edit(hp, [edit_pos_x edit_pos_y-120 80 30], tooltip, @valid_positive_number);
    
    gui_lib.uicontrol_text(hp, 'Relative change:', [title_pos_x title_pos_y-160 150 30]);
    val_options = {'no', 'difference', 'percentage'};
        
    tooltip = '<html><b>Relative change</b><br><i>Select to desplay relative changes: <li> difference <li> percentage</html>';
    relative = uicontrol(hp,'Style','popupmenu',...
        'String', val_options, 'Position',[edit_pos_x-40 edit_pos_y-160 120 25], 'TooltipString', tooltip);%samples pre-event

    gui_lib.uicontrol_text(hp, 'Baseline:', [title_pos_x title_pos_y-240 150 30]);
    tooltip = '<html><b>Beseline</b><br><i>Select number of ms before the first event for defining the baseline.';
    baseline = gui_lib.uicontrol_edit(hp, [edit_pos_x edit_pos_y-240 80 30], tooltip, @valid_positive_number);
    
    val_options = {'no', 'SE', 'CI'};

    gui_lib.uicontrol_text(hp, 'Scattering:', [title_pos_x title_pos_y-280 150 30]);
    tooltip = '<html><b>Scattering</b><br><i>Select scattering approach around the mean. <li> SE - Starndard Error <li> CI - Confidence Interval ';
    scattering = uicontrol(hp,'Style','popupmenu',...
        'String', val_options, 'Position',[edit_pos_x edit_pos_y-280 80 30], 'TooltipString', tooltip);%samples pre-event


    
 
    tooltip = '<html><b>Method</b><br><i>Select method: <li> Arbitrary units <li> Z score ';
    gui_lib.uicontrol_text(hp, 'Method:', [title_pos_x title_pos_y-200 150 30]);
    val_options = {'arbitrary units'; 'Z score'};
    Method =  uicontrol(hp,'Style','popupmenu',...
        'String', val_options, 'Position',[edit_pos_x-40 edit_pos_y-200 120 25], 'TooltipString', tooltip);

    uistack(to, 'up', 3);
    uicontrol(BinsNumber);

    data.configuration.BinsNumber         = BinsNumber;
    data.configuration.PreEventNumber     = PreEventNumber;
    data.configuration.Method             = Method;

    data.configuration.from               = from;
    data.configuration.to                 = to;
    data.configuration.baseline           = baseline;
    data.configuration.relative           = relative;
    data.configuration.scattering         = scattering;

end

function handler = view_user_vars(figure, conds_data, vars)
    handler = [];

    buildin_var = {'pupil_x', 'pupil_y', 'pupil_size', 'trial_id', 'event_Trial_Offset', 'timestamps', 'msgs'};
    vars        = vars(~ismember(vars, buildin_var));
    
    hp_vars     = uipanel(figure, 'Units', 'Pixels', 'BackgroundColor', [0.1 0.32 0.46], 'BorderWidth', 0, 'Position', [285 374 270 270]);   
    gui_lib.uicontrol_title('USER''S VARIABLES', 345, 95);

    hp_events = uipanel('Units','Pixels', 'BackgroundColor',[0.1 0.32 0.46], 'BorderWidth', 0, 'Position', [595 374 270 270]);
    gui_lib.uicontrol_title('USER''S EVENTS', 670, 95);

    if (isempty(vars))
       return
    end
    
    events  = vars(~cellfun(@isempty, strfind(vars,'event_')));
    vars    = vars(cellfun(@isempty, strfind(vars,'event_')));
    events  = cellfun(@(x) x(7:end), events, 'UniformOutput', false); 
    text_pos_x   = 10;
    text_pos_y   = 240;
    max_var2show = 8;
    max_buts     = 8;

    handler = gui_lib.paging(hp_vars, handler, vars,  [text_pos_x text_pos_y 250 30], max_var2show, max_buts);
    handler = gui_lib.paging(hp_events, handler, events,  [text_pos_x text_pos_y 250 30], max_var2show, max_buts, 'event');

    % presents variable's values
    gui_lib.uicontrol_title('VARIABLE VALUES', 380, 435);
    
    width = 812;
    hight = 150;
    hp_table = uipanel('Units','Pixels', 'BackgroundColor',[0.1 0.32 0.46],...
                  'Position',[44 150 width  hight]);
    y     = 147-hight;
    x     = 1;
    rows  = 1:size(conds_data,1);
    
    uitable(hp_table, 'Data', conds_data, 'position', [x y width-1 hight-1],...
                 'ColumnName', char(vars),'RowName', rows, 'ColumnWidth',{100}); 
end


function view_analyze_results(analyzed_data, comp_names, data_table, src, cond_mat_cleaned, cond_events, figure, fig, log, log_a)
    width = 822;
    hight = 150;
    
    hpt   = uipanel(figure, 'Units','Pixels', 'BackgroundColor',[0.1 0.32 0.46], 'Position',[38 90 width  150]);
    y     = 147-hight;
    x     = 1;

    header = [cellstr('Outliers') cellstr('Valid_Trials') fieldnames(analyzed_data.(char(comp_names(1))).events)'];
    header = cellfun(@(x) iff(strcmp(x(1:6),'event_'),  x(7:end) ,x), header, 'UniformOutput', false);
    header = strrep(header, '_', ' ');
    comp_names_fixed = cellfun(@(x) x(3:end), comp_names, 'UniformOutput', false);
    comp_names = strrep(strrep(comp_names_fixed, '_x_', ' & '),'_',' ');
 
    res_table = uitable(hpt, 'Data', data_table, 'position', [x y width-1 hight-1],...
                      'ColumnName', header','RowName', comp_names, 'ColumnWidth',{100}); 
    
    gui_lib.uicontrol_button(figure, [32 40 125 30], 'Exit', @exit);
    tooltip = '<html><i>Select again the conditions for the analysis';

    gui_lib.uicontrol_button(figure, [162 40 155 30],'Change conditions', {@(~,~)change_comps(src, cond_mat_cleaned, cond_events)}, tooltip); 
    tooltip = '<html><i>Select folder with data files and <br> run the same parameters on them';
    gui_lib.uicontrol_button(figure, [322 40 125 30],'Group analysis', {@(~,~)do_group(src, fig, log, log_a, res_table)}, tooltip); 

    tooltip = '<html><i>Save the figure of the participant as Matlab figure (fig file)';
    gui_lib.uicontrol_button(figure, [354 281 155 30], 'Save figure', @(~,~)output.save_figure(fig, comp_names_fixed), tooltip);
%     tooltip = '<html><i>Save the figure of the participant as image (png file)';
%     gui_lib.uicontrol_button(figure, [534 281 155 30], 'Save as image', @(~,~)output.(fig, comp_names_fixed), tooltip);
    tooltip = '<html><i>Save the data of the participant (mat file)';
    gui_lib.uicontrol_button(figure, [714 281 155 30], 'Save data', @(~,~)save_participant_output(src, fig, log, log_a), tooltip);   

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  get/set functions                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function configuration = Get_processing_vals(configuration)

    min_trials = get(configuration.min_trials,'String');
    configuration.min_trials_val = str2double(min_trials);

    ZNumber = get(configuration.ZNumber,'String');
    configuration.ZOutliers_val = str2double(ZNumber);

    ZeroshNumber = get(configuration.ZeroshNumber,'String');
    configuration.ZeroshNumber_val = str2double(ZeroshNumber);
end

function save_participant_output(src, fig, log, log_a)
    data = guidata(src);
    do_analyze(data, fig, src, log, log_a);
    output.save_output(src)
end

function configuration = get_analyze_vals(configuration)

    configuration.BinsNumber_val          = str2double(get(configuration.BinsNumber, 'String'));
    if(isnan(configuration.BinsNumber_val) || configuration.BinsNumber_val<=0)
        configuration.BinsNumber_val = 1;
    end    
    
    configuration.relative_val            = gui_lib.get_popupmenu_val(configuration.relative);
    
    configuration.Method_val              = gui_lib.get_popupmenu_val(configuration.Method);    

    configuration.scattering_val          = gui_lib.get_popupmenu_val(configuration.scattering);    

    configuration.from_val                = gui_lib.get_popupmenu_val(configuration.from);
    configuration.to_val                  = gui_lib.get_popupmenu_val(configuration.to);
    configuration.comp_names_val          = configuration.comp_names;
    
    configuration.PreEventNumber_val      = str2double(get(configuration.PreEventNumber, 'String'));
    if isnan(configuration.PreEventNumber_val)
        configuration.PreEventNumber_val = 0;
    end

    configuration.baseline_val  = str2double(get(configuration.baseline, 'String'));

    if isnan(configuration.baseline_val)
        configuration.baseline_val = 0;
    end
end

function var_handler   = save_vars(handler, levels)
    events = [];
    conditions = [];
    var_names     = fieldnames(handler);
    for i=1:size(var_names, 1)
        var = handler.(char(var_names(i)));
        val = gui_lib.get_checkbox_val(var);
        if(val>0)
            if(isempty(strfind(char(var_names(i)), 'event_')))
                cond_levels = levels(:, i);
                cond_levels = cond_levels(~strcmp(cond_levels, ''));
                conditions.(char(var_names(i))) = cond_levels; 
            else
                events = [events var_names(i)];
            end
        end
    end
    var_handler.events     = events;
    var_handler.conditions = conditions;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   utiles functions                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = iff(condition, trueResult, falseResult)
    if condition
        result = trueResult;
    else
        result = falseResult;
    end
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

function valid_positive_number(src, ~)
    num = str2num(get(src, 'string'));
    if  isempty(num) || num<0
        set(src,'string',0)
    end
end


function valid_2plus_number(src, ~)
    num = str2num(get(src, 'string'));
    if  isempty(num) || num<2
        set(src,'string',2)
    end
end
