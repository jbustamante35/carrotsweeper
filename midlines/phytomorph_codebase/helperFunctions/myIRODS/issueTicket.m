function [fileName iticket] = issueTicket(fileName,uses)
    if nargin == 1
        uses = 1;
    end
    cmd = ['iticket create read "' fileName '"'];
    [o,r] = system(cmd);
    fidx = strfind(r,':');
    iticket = r(fidx+1:end-1);
    cmd = ['iticket mod ' iticket ' uses ' uses];
    [o,r] = system(cmd);
    fileName = [fileName '#' iticket '#'];
end