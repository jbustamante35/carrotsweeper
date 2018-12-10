function [FileList] = issueBulkTicket(FileList)
    ticketMultiplier = 2;
    if ~isempty(FileList)
        
        
        % get the paths to the file names
        for e = 1:numel(FileList)
            [pth{e},nm,ext] = fileparts(FileList{e});
        end
        UQ = unique(pth);
        
        
        % find the directories the files belong to
        for u = 1:numel(UQ)
            ticketCNT(u) = 0;
            for e = 1:numel(pth)
                if strcmp(pth{e},UQ{u})
                    directoryU(e) = u;
                    ticketCNT(u) = ticketCNT(u) + 1;
                end
            end
        end
        
        
        for u = 1:numel(UQ)
            fidx = find(directoryU == u);
            [junkName iticket] = issueTicket(UQ{u},ticketMultiplier*ticketCNT(u));
            for e = 1:numel(fidx)
                FileList{fidx(e)} = [FileList{fidx(e)} '#' iticket '#'];
            end
        end
    end
end