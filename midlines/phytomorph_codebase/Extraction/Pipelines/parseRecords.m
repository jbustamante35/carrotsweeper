function [r] = parseRecords(text)
    if strcmp(text(1),char(10))
        text(1) = [];
    end
    rs = ['------------------------------------------------------------' char(10)];
    %rs = ['------------------------------------------------------------'];
    text = [rs text];
    fidx = strfind(text,rs);
    cnt = 1;
    for e = 1:(numel(fidx)-1)
        try
            record = text((fidx(e)+numel(rs)):(fidx(e+1)-1));
            record = [char(10) record];
            gidx = strfind(record,'=');
            lidx = strfind(record,char(10));
            for f = 1:(numel(lidx)-1)
                key = record((lidx(f)+1):(gidx(f)-2));
                value = record(gidx(f)+2:(lidx(f+1)-1));
                r(cnt).(key) = value;
            end
            cnt = cnt + 1;
        catch
        end
    end
end