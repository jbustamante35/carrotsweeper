function [result] = ils(path)    
    [status,result] = system('ils');
    %{
    if nargin == 0;
        [status path] = system('ipwd');
    end
    %[status, path] = system(path);
    %fidx_end = strfind(result,sprintf('\n'));
    
    cmd = ['iquest "SELECT DATA_ID, DATA_NAME WHERE COLL_NAME = ''' path(1:end-1) ''' "'];            
    [status query] = system(cmd);            
    fidx_end = strfind(query,sprintf('\n'));
    record = 1;
    for i = 1:(size(fidx_end,2)-1)
        line{i} = query(fidx_end(i)+1:fidx_end(i+1)-1);
        if strcmp(line{i},'------------------------------------------------------------')
            record = record + 1;
        end
        fidx = strfind(line{i},'=');
        if ~isempty(fidx)
            fld = strtrim(line{i}(1:fidx-1));
            S(record).(fld) = strtrim([line{i}(fidx+1:end)]);
            %S(record).value = strtrim([line{i}(fidx+1:end)]);
        end
    end
    %}
    
    %{
    for i = 1:(size(fidx_end,2)-1)
        line{i} = result(fidx_end(i)+1:fidx_end(i+1)-1);
        tmp = strfind(line{i},'C- ');
        if ~isempty(tmp)
            S(i).name = [path(1:end-1) '/' line{i}(4:end)];
            S(i).isdir = logical(1);            
            % get the unique ID 
            cmd = ['iquest "SELECT COLL_ID WHERE DATA_NAME = ''' line{i} ''' AND COLL_NAME = ''' path(1:end-1) ''' "'];            
            [status query] = system(cmd);            
            % store the unique ID
            fidx0 = strfind(query,'DATA_ID');
            fidx1 = strfind(query,sprintf('\n'));
            nidx = find(fidx1 > fidx0);            
            S(i).DATA_ID = strtrim(query(fidx0+9:fidx1(nidx(1))-1));
        else
            line{i} = strtrim(line{i});
            S(i).name = [path(1:end-1) '/' line{i}(3:end)];
            S(i).isdir = logical(0);
            % get the unique ID 
            cmd = ['iquest "SELECT DATA_ID WHERE DATA_NAME = ''' line{i} ''' AND COLL_NAME = ''' path(1:end-1) ''' "'];            
            [status query] = system(cmd);            
            % store the unique ID
            fidx0 = strfind(query,'DATA_ID');
            fidx1 = strfind(query,sprintf('\n'));
            nidx = find(fidx1 > fidx0);            
            S(i).DATA_ID = strtrim(query(fidx0+9:fidx1(nidx(1))-1));
        end
        i
    end
    %}
end