function print_log(log_str, log)
    if (~exist('log', 'var') || strcmp(log, ''))
        fprintf('%s\n', log_str);
        return;
    end
    
    set(log, 'String', [' ' log_str]);
    drawnow
    if strcmp(log_str, '')
        return;
    end
    
    full_path = mfilename('fullpath');
    path      = full_path(1:end-length(mfilename));
    log_path  = strcat(path, filesep, 'logs');
    if ~exist(log_path, 'dir')
        mkdir(log_path);
    end
    
    fid = fopen(strcat(log_path, filesep, date(), '.txt'), 'a+') ;
    time = datestr(str2double(char( num2str(now, '%.11f'))), 'yyyy-mm-dd HH:MM:SS.FFF');
    fprintf(fid, '%s%s%s\n', time, ' | ', log_str);
    fclose(fid) ;
end
