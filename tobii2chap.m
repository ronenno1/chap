[file_name, path]  = uigetfile('*.csv', 'Select csv file ');

if ~file_name
    return
end

relevant_fields = {'PupilDiameterLeft'; 'PupilDiameterRight'; 'RecordingTimestamp'};

[params_file, params_path]  = uigetfile('*.csv', 'Select params file ');
% The file should includes these headers: trial_onset, times, vars, events, pupil_left, pupil_right
% 
% For example:    trial_onset	times               vars	events          pupil_left          pupil_right
%                 TrialId       RecordingTimestamp	ACC     CurrentObject	PupilDiameterLeft   PupilDiameterRight
%                                                   block		
%                                                   TASK

if ~params_file
    return
end

full_parame_file = [params_path, params_file];
tic
params    = readtable(full_parame_file);

full_csv_name = [path, file_name];
[~, file_name, ~] = fileparts(full_csv_name);
full_name = [path, file_name];


trial_times = params.times{1, :};

trial_onset = params.trial_onset{1, :};


vars = params.vars;
vars =  vars(~cellfun('isempty', vars));

pupil_left  = params.pupil_left{1, :};
pupil_right = params.pupil_right{1, :};

tobii_file = readtable([full_name '.csv']);

events_input = tobii_file.(params.events{1});

[events, pos] = unique(events_input);
[~, idx] = sort(pos);
events = events(idx);
events =  events(~cellfun('isempty', events));


[lo, trial_onset_ids] = unique(tobii_file.(trial_onset));
trial_onset_ids = sort(trial_onset_ids);
trial_onset_ids = trial_onset_ids+1;
% trial_onset_ids(1) = trial_onset_ids(1)+1;
all_times = tobii_file.(char(trial_times));

all_times_sort = sort(all_times);
if sum(all_times_sort ~= all_times)>0
    find((all_times_sort ~= all_times)>0, 1, 'first')
    disp('There is a problem with the times');
    return
end

% onsets
disp('Parsing trial onsets...');

trial_onsets = all_times(trial_onset_ids);
events_data.TRIALID = trial_onsets;



% offsets
disp('Parsing trial offsets...');

ms = round(mean(diff(tobii_file.(trial_times)(1:50))));
trial_offsets = zeros(size(trial_onsets));
for trial = 1:length(trial_onset_ids)-1
    trial_offsets(trial) = trial_onsets(trial+1)-ms;
end
trial_offsets(end) = tobii_file.(trial_times)(end);


% vars
disp('Parsing variables...');

for var = 1:length(vars)
    var_data.(vars{var}) = tobii_file.(vars{var})(trial_onset_ids);
end


disp('Parsing events...');

events_struct.Event = {};
events_struct.RecordingTimestamp = [];

for trial = 1:length(trial_onset_ids)
    events_struct.RecordingTimestamp(end+1, :) = trial_onsets(trial);
    events_struct.Event(end+1, :) = {['TRIALID ' num2str(trial)]};
    for event = 1:length(events)
        event_time = trial_onset_ids(trial) + find(~cellfun(@isempty, strfind(events_input(trial_onset_ids(trial):end), events{event})), 1);
        events_struct.Event(end+1, :) = {['!E TRIAL_EVENT_VAR ' events{event}]};
        events_struct.RecordingTimestamp(end+1, :) = tobii_file.(trial_times)(event_time);        
    end
    
    for var = 1:length(vars)
        events_struct.RecordingTimestamp(end+1, :) = trial_offsets(trial);    
        val = var_data.(vars{var})(trial);
        if ~iscell(val)
            val = num2str(val);
        end
        events_struct.Event(end+1, :) = strcat('!V TRIAL_VAR', {' '}, vars(var), {' '}, val);
    end


    events_struct.RecordingTimestamp(end+1) = trial_offsets(trial);
    events_struct.Event(end+1, :) = {'TRIAL_END'};
end
events_table = struct2table(events_struct);
disp('Writing files...');

writetable(events_table, [full_name, '_events.csv']);

% prepare and save tbi file 
headers  = tobii_file.Properties.VariableNames;
header_ids = zeros(3, 1);
[~, header_ids(1)] = find(strcmp(headers, trial_times));
headers(header_ids(1)) = {'RecordingTimestamp'};


[~, header_ids(2)] = find(strcmp(headers, pupil_left));
headers(header_ids(2)) = {'PupilDiameterLeft'};

[~, header_ids(3)] = find(strcmp(headers, pupil_right));
headers(header_ids(3)) = {'PupilDiameterRight'};


tobii_file.Properties.VariableNames = headers;

tobii_file.PupilDiameterLeft(tobii_file.PupilDiameterLeft==-1) = 0;
tobii_file.PupilDiameterRight(tobii_file.PupilDiameterRight==-1) = 0;

headers_ids2remove = ones(1, length(headers));
headers_ids2remove(header_ids) = 0;
tobii_file = removevars(tobii_file, headers(find(headers_ids2remove)));

writetable(tobii_file, [full_name, '_new.csv']);
copyfile([full_name, '_new.csv'], [full_name, '.tbi']);    
delete([full_name, '_new.csv']);
a = toc;
disp(['Done (elapsed time is ', num2str(toc), ' seconds)!']);
