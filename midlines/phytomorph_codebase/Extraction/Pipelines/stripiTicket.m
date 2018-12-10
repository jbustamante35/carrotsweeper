function [rawName,ticket] = stripiTicket(fileName)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % strip ticket data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['***********************************\n']);
    fprintf(['Start ticket strip:' fileName '\n']);
    fprintf(['***********************************\n']);
    if ~isempty(fileName)
        fidx = strfind(fileName,'#');
        if ~isempty(fidx)
            rawName = fileName(1:(fidx(1)-1));
            ticket = fileName((fidx(1)+1):(fidx(2)-1));
        else
            rawName = fileName;
            ticket = '';
        end
    else
        rawName = fileName;
        ticket = '';
    end
    fprintf(['***********************************\n']);
    fprintf(['End ticket strip:' fileName '\n']);
    fprintf(['Raw Name-->' rawName '\n']);
    fprintf(['Ticket Data-->' ticket '\n']);
    fprintf(['***********************************\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % strip ticket data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end