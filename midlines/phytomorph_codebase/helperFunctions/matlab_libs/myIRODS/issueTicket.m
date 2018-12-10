function [fileName iticket] = issueTicket(fileName,uses,type)
    if nargin == 1
        uses = 1;
    end
    if nargin <= 2
        type = 'read';
    end
    cmd = ['iticket create ' type ' ' '"' fileName '"'];
    [o,r] = system(cmd);
    fidx = strfind(r,':');
    iticket = r(fidx+1:end-1);
    if uses ~= Inf
        cmd = ['iticket mod ' iticket ' uses ' num2str(uses)];
    end
    %[o,r] = system(cmd);
    if strcmp(type,'write')
        cmd = ['iticket mod ' iticket ' write-file ' num2str(uses)];
        cmd = ['iticket mod ' iticket ' write-file ' num2str(0)];
        [o,r] = system(cmd);
    end
    expireDate = datestr(addtodate(now,11,'day'),'YYYY-mm-dd.hh:MM:ss');
    cmd = ['iticket mod ' iticket ' expire ' expireDate];
    [o,r] = system(cmd);
    fileName = [fileName '#' iticket '#'];
end