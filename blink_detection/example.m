load('trials_data.mat');                        % load test samples data
sampling_rate          = 500;                   % set the sampling rate
sampling_interval      = 1000/sampling_rate;	% compute the sampling time interval in milliseconds.
output_folder          = 'output';              % set the output directory
show_and_save_figures  = true;                  % control whether to show & save result plots

if ~exist(output_folder, 'dir')	% if output directory does not exist - create it
    mkdir(output_folder);
end;

fid = fopen(strcat(output_folder, filesep, 'blinks.csv'), 'w') ;	% open results file
fig = figure(1);                                                    % create new figure for plotting

for trial_id = 1:size(trials_data, 1)
    trial_data = trials_data(trial_id, :);                          % extract trial data
    trial_data(isnan(trial_data)) = [];                             % remove irrelevant NaNs at the end of the trial

    blinks_data_positions     = based_noise_blinks_detection(trial_data', sampling_rate);	% get blink positions using the noise-based approach
    blinks_data_positions_str = strjoin(arrayfun(@(x) num2str(x), blinks_data_positions' ,'UniformOutput', false), ',');	% format results as a string for printing
    fprintf(fid, blinks_data_positions_str) ;	% print blink positions
    fprintf(fid, '\n') ;    

    %% plot and save results
    if(~show_and_save_figures)	% if there is no need for displaying and saving the result plots, move to the next trial
        continue;
    end;
    cla(fig);
    x_axis = linspace(0 ,sampling_interval*size(trial_data, 2), size(trial_data, 2));
    trial_data(trial_data==0) = nan;
    plot(x_axis, trial_data);
    hold on
    plot(blinks_data_positions, trial_data(blinks_data_positions/sampling_interval), 'ro', 'LineWidth', 3);
    xlabel('Time [ms]', 'FontWeight','bold');
    ylabel('Pupil Size [arbitrary units]', 'FontWeight','bold');
    set(gca, 'FontWeight', 'bold');
    title(strcat('Trial: ', num2str(trial_id))); 
    saveas(fig, strcat(output_folder, filesep, num2str(trial_id)));	                  % save plot to output directory as a figure file
    print(fig, '-dpng', '-r300', strcat(output_folder, filesep, num2str(trial_id)));  % save plot to output directory as a png file
end;
close(fig);