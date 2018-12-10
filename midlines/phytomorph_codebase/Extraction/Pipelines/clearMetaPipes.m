function [] = clearMetaPipes(user,plantType,tissueType)
    %%
    tic
    CMD = ['iquest --no-page "SELECT DATA_NAME,COLL_NAME,META_DATA_ATTR_VALUE,META_DATA_ATTR_NAME,META_DATA_ATTR_UNITS,META_DATA_MODIFY_TIME ' ... 
           'WHERE COLL_NAME LIKE ''/iplant/home/' user ...
           '/' plantType 'Data/' tissueType 'Data%''"'];
    [o,r] = system(CMD);
    fprintf(['Done running SQL string:' num2str(toc) '\n']);
    %%
    [r] = parseRecords(r);
    %%
    tic
    rm = [];
    for e = 1:numel(r)
        if strcmp(r(e).META_DATA_ATTR_NAME,'ipc_UUID')
            rm = [rm e];
        end
    end
    r(rm) = [];
    fprintf(['Done parsing resulting SQL statement:' num2str(toc) '\n']);
    %%
    tic
    string = '';
    tmpShell = [tempname '.r'];
    fileID = fopen(tmpShell,'w');
    for e = 1:numel(r)
        %string = [string [r(e).COLL_NAME filesep r(e).DATA_NAME] char(10) r(e).META_DATA_ATTR_NAME char(10) r(e).META_DATA_ATTR_VALUE char(10)];
        if isempty(r(e).META_DATA_ATTR_UNITS)
            eS = char(10);
        else
            eS = [' "' r(e).META_DATA_ATTR_UNITS '"' char(10)];
        end
        string = [' rm -d "' r(e).COLL_NAME filesep r(e).DATA_NAME '" "' r(e).META_DATA_ATTR_NAME '" "' r(e).META_DATA_ATTR_VALUE '"' eS];
        fprintf(fileID,'%s',string);
    end
    
    % write to shell script
   
    fclose(fileID);
    fprintf(['Done generating removal string:' num2str(toc) '\n']);
    tic;
    cmd = ['parallel --gnu -u -j 20 -n 1 -a ' tmpShell ' myImeta.sh '];
    tic
    [o,res] = system(cmd);
    fprintf(['Done running removal string:' num2str(toc) '\n']);
end

%{
 clearMetaPipes('gxe','maize','cob')
%}