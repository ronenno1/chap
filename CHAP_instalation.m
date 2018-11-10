function CHAP_instalation()
clear mex

targetdirectory=fileparts(mfilename('fullpath'));
cd (targetdirectory);

% Check OS
IsWin = ~isempty(strfind(computer, 'PCWIN')) || ~isempty(strfind(computer, '-w64-mingw32'));
IsOSX = ~isempty(strfind(computer, 'MAC')) || ~isempty(strfind(computer, 'apple-darwin'));
IsLinux = strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64') || ~isempty(strfind(computer, 'linux-gnu'));



% Does SAVEPATH work?
if exist('savepath')
   err=savepath;
else
   err=path2rc;
end

if err
    try
        % If this works then we're likely on Matlab:
        p=fullfile(matlabroot,'toolbox','local','pathdef.m');
        fprintf(['Sorry, SAVEPATH failed. Probably the pathdef.m file lacks write permission. \n'...
                 'Please ask a user with administrator privileges to enable \n'...
                 'write by everyone for the file:\n\n''%s''\n\n'],p);
    catch
        % Probably on Octave:
        fprintf(['Sorry, SAVEPATH failed. Probably your ~/.octaverc file lacks write permission. \n'...
                 'Please ask a user with administrator privileges to enable \n'...
                 'write by everyone for that file.\n\n']);
    end
    
    fprintf(['Once "savepath" works (no error message), run ' mfilename ' again.\n']);
    fprintf('Alternatively you can choose to continue with installation, but then you will have\n');
    fprintf('to resolve this permission isssue later and add the path to the CHAP manually.\n\n');
    answer=input('Do you want to continue the installation despite the failure of SAVEPATH (yes or no)? ','s');
    if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
        fprintf('\n\n');
        error('SAVEPATH failed. Please get an administrator to allow everyone to write pathdef.m.');
    end
end

searchpattern = [pwd '[' filesep pathsep ']'];
searchpattern2 = pwd;


% Remove "CHAP" from path:
while any(regexp(path, searchpattern))
    fprintf('Your old CHAP appears in the MATLAB/OCTAVE path:\n');
    paths=regexp(path,['[^' pathsep ']+'],'match');
    answer=input('Before you decide to delete the paths, do you want to see them (yes or no)? ','s');
    if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
        fprintf('You didn''t say "yes", so I''m taking it as no.\n');
    else
        for p=paths
            s=char(p);
            if any(regexp(s,searchpattern2))
                fprintf('%s\n',s);
            end
        end
    end
    answer=input('Shall I delete all those instances from the MATLAB/OCTAVE path (yes or no)? ','s');
    if ~strcmpi(answer,'yes') && ~strcmpi(answer,'y')
        fprintf('You didn''t say yes, so I cannot proceed.\n');
        fprintf('Please use the MATLAB "File:Set Path" command or its Octave equivalent to remove all instances of "CHAP" from the path.\n');
        error('Please remove CHAP from MATLAB/OCTAVE path.');
    end
    for p=paths
        s=char(p);
        if any(regexp(s,searchpattern2))
            rmpath(s);
        end
    end
    if exist('savepath') %#ok<EXIST>
       savepath;
    else
       path2rc;
    end

    fprintf('Success.\n\n');
end

% if ~strcmpi(targetdirectory, pwd)
%     targetdirectory=fileparts(mfilename('fullpath'));
% end


% Add CHAP to MATLAB / OCTAVE path
fprintf('Now adding the new CHAP folder (and all its subfolders) to your MATLAB / OCTAVE path.\n');

p=fullfile(targetdirectory);
pp=genpath(p);
pathElements = textscan(pp, '%s', 'delimiter', pathsep);
pathElements = pathElements{1}.';

isNotMatching    = cellfun(@isempty,strfind(pathElements, '.git'));
pathElements     = pathElements(isNotMatching);
pathElements     = [pathElements; repmat({pathsep},1,length(pathElements))];
newPathList      = [pathElements{:}];
newPathList(end) = [];

addpath(newPathList);


if exist('savepath')
   err=savepath;
else
   err=path2rc;
end

if err
    fprintf('SAVEPATH failed. CHAP is now already installed and configured for use on your Computer,\n');
    fprintf('but i could not save the updated MATLAB / OCTAVE path, probably due to insufficient permissions.\n');
    fprintf('You will either need to fix this manually via use of the path-browser (Menu: File -> Set Path),\n');
    fprintf('or by manual invocation of the savepath command (See help savepath). The third option is, of course,\n');
    fprintf('to add the path to the CHAP folder and all of its subfolders whenever you restart MATLAB / OCTAVE.\n\n\n');
else 
    fprintf('Success.\n\n');
end

if ismac
    zip_file = [pwd, filesep, 'edf2mat', filesep, 'edf-converter-master.zip'];
    zip_dest = [pwd, filesep, 'edf2mat'];
    folder_1 = [pwd, filesep, 'edf2mat', filesep, 'edf-converter-master', filesep, '@Edf2Mat'];
    folder_2 = [pwd, filesep, 'edf2mat', filesep, 'edf-converter-master', filesep, 'edfmex'];
    dest_1   = [pwd, filesep, '@Edf2Mat'];
    dest_2   = [pwd, filesep, 'edfmex'];
    try
        unzip(zip_file, zip_dest); 
        copyfile(folder_1, dest_1);
        copyfile(folder_2, dest_2);
    catch
        fprintf('Something works wrong :(\nPlease unzip: %s\n', zip_file);
        fprintf('Then copy: %s and %s to %s\n', folder_1, folder_2, pwd);
    end
end
