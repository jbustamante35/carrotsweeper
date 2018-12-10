function [S] = iquery2struct(query)
    % frame the query
    recSep = '------------------------------------------------------------';
    query = [recSep sprintf('\n') query];
    fidx_record = strfind(query,'------------------------------------------------------------');
    % cycle through the frames
    for i = 1:(size(fidx_record,2)-1)
        rec{i} = query(fidx_record(i)+numel(recSep):fidx_record(i+1)-1);
        equal_fidx = strfind(rec{i},'=');
        newline_fidx = strfind(rec{i},sprintf('\n'));        
        for j = 1:size(equal_fidx,2)
            field_name = strtrim(rec{i}(newline_fidx(j)+1:equal_fidx(j)-1));
            value = strtrim(rec{i}(equal_fidx(j)+1:newline_fidx(j+1)-1));
            S(i).(field_name) = value;
        end
    end
end