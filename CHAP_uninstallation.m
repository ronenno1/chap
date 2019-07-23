function CHAP_uninstallation()
clear mex

targetdirectory=fileparts(mfilename('fullpath'));
cd (targetdirectory);


% Does SAVEPATH work?
if exist('savepath')
   err=savepath;
else
   err=path2rc;
end

if err
    p=fullfile(matlabroot,'toolbox','local','pathdef.m');
    fprintf(['Sorry, SAVEPATH failed. Probably the pathdef.m file lacks write permission. \n'...
             'Please ask a user with administrator privileges to enable \n'...
             'write by everyone for the file:\n\n''%s''\n\n'],p);
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

if ~any(regexp(path, searchpattern))
    fprintf('CHAP is not installed\n');
    return 
end
% Remove "CHAP" from path:
while any(regexp(path, searchpattern))
    paths=regexp(path,['[^' pathsep ']+'],'match');
    for p=paths
        s=char(p);
        if any(regexp(s,searchpattern2))
            rmpath(s);
            fprintf('Remove from MATLAB path: %s \n',s);
            
        end
    end
    if exist('savepath') %#ok<EXIST>
       savepath;
    else
       path2rc;
    end
    fprintf('CHAP successfully removed.\n\n');
end