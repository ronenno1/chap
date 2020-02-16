% In this example we demonstrate how to add a variable to an existing chp file.
% Specifically, in the current example we demonstrate who to use a list of 
% participants that should be excluded from the pupillometry analysis 
% (for example, because they failed to perform the task correctly).
% The list should include one field with the header 'id'.
% The script will add the value ‘invalid’ to participants who appear in the list, 
% and the value ‘valid’ to all other participants.


%% getting a list of all the chp files
chp_files = dir(['*chp']);
chp_files = {chp_files.name}';

%% reading the list - this list is relevant only for the specific example
list = readtable('list.csv');

%% run across all the chp files
for id = 1:length(chp_files)
    %% get the file name and load ot
    [~, sub_id, ~] = fileparts(chp_files{id});
    disp(['Fixing ' sub_id '...']);
    sub = load(chp_files{id}, '-mat');

    %% here is the logic of the specific example
    if ~isempty(find(strcmp(list.id, sub_id)))      % check if the subject is included in the list
        validity = 'invalid';                       % the subject is invalid
    else
        validity = 'valid';                         % the subject is valid
    end
    
    %% adding the new variable 
    for i = 1:size(sub.data.total_var_data_table, 1)    
        sub.data.total_var_data_table.validity{i} = validity;
    end

    %% saving the data including the new variable
    data = sub.data;
    save([sub_id '.chp'], 'data');
end