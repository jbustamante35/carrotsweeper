function [S] = iquest(info,cond,rest)
    cmd = ['iquest "SELECT '];
    % build selection
    for i = 1:size(info,2)
        cmd = [cmd info{i} ','];
    end
    cmd(end) = [];
    cmd = [cmd ' WHERE '];
    % build conditions
    for i = 1:size(cond,2)
        cmd = [cmd cond{i} '=''' rest{i} ''' AND '];
    end
    cmd(end-3:end) = [];
    cmd = [cmd '"'];
    
    [status query] = system(cmd);            

end