function [CMD] = uncLog(preFix,fileList,cmd,algo,ver,stage,value,execute)


    % if no stages are entered - then fill all
    % else fill N < all
    if isempty(stage)
        N = 6;
    else
        N = numel(stage);
    end
    switch cmd
        case 'add'
            %{
            for s = 1:N
                % if no stages are given
                % fill all
                
                stage{s} = num2str(s);

                value{s} = '0';
            end
            %}
        case 'rmw'
            for s = 1:N
                stage{s} = num2str(s);
                value{s} = '%';
            end
    end
    
    if numel(fileList) > 1
        tic
        tmpShell = [tempname '.r'];
        fileID = fopen(tmpShell,'w');
        
        for e = 1:numel(fileList)
            tmp = generateForFile(preFix,fileList{e},cmd,algo,ver,stage,value,'stream');
            for s = 1:numel(tmp)
                fprintf(fileID,'%s',[' ' tmp{s} char(10)]);
            end
        end
        % write to shell script
        fclose(fileID);
        
        
        fprintf(['Done generating removal string:' num2str(toc) '\n']);
        tic;
        CMD = ['parallel --gnu -u -j 20 -n 1 -a ' tmpShell ' myImeta.sh '];
        tic
        [o,res] = system(CMD);
        fprintf(['Done running string:' num2str(toc) '\n']);
        
        %{
        
        
        %CMDs = ['imeta << EOF' char(10)];
        tempLate1 = 'msiString2KeyValPair(#KVP#,*var);';
        if strcmp(cmd,'add')
            tempLate2 = 'msiAssociateKeyValuePairsToObj(*var, "#filename#", "-d");';
        elseif strcmp(cmd,'rmw')
            tempLate2 = 'msiRemoveKeyValuePairsFromObj(*var, "#filename#", "-d");';
        end
        
        
        tmp = generateForFile(preFix,[],cmd,algo,ver,stage,value,'stream');
        wholeKVP = '''';
        for i = 1:numel(tmp)
            wholeKVP = [wholeKVP tmp{i} '%'];
        end
        wholeKVP(end) = '''';
        tempLate1 = strrep(tempLate1,'#KVP#',wholeKVP);
        
        
        CMDs = ['tempFunction { ' char(10) tempLate1 char(10)];
        for e = 1:numel(fileList)
            CMDs = [CMDs strrep(tempLate2,'#filename#',fileList{e}) char(10)];
        end
        CMDs = [CMDs char(10) '}' char(10) 'INPUT null' char(10) 'OUTPUT ruleExecOut'];
        CMD{1} = CMDs;
        %}
    else
        CMD = generateForFile(preFix,fileList{1},cmd,algo,ver,stage,value,'single');
    end
    
    
    if execute
        for e = 1:numel(CMD)
            if numel(CMD) > 1
                %{
                tm = clock;
                fprintf(['Start generating shell script:\n']);
                % create temp shell script
                tmpShell = [tempname '.r'];
                fileID = fopen(tmpShell,'w');
                % write to shell script
                fprintf(fileID,'%s\n',CMD{e});
                fclose(fileID);
                fprintf(['Done generating shell script:' num2str(etime(clock,tm)) '\n']);
                tm = clock;
                fprintf(['Start running iMeta:' cmd '\n']);
                cmd = ['irule -F ' tmpShell];
                [r,o] = system(cmd);
                delete(tmpShell);
                fprintf(['Done running iMeta:' cmd ':'  num2str(etime(clock,tm)) '\n']);
                %}
            else
                tm = clock;
                cmd = CMD{e};
                fprintf(['Start running iMeta:' cmd '\n']);
                [r,o] = system(cmd);
                fprintf(['Done running iMeta:' cmd ':'  num2str(etime(clock,tm)) '\n']);
            end
        end
    end
    
    
end

function [CMD] = generateForFile(preFix,file,cmd,algo,ver,key,value,streamType)
    CMD = {};
    if strcmp(cmd,'add') || strcmp(cmd,'rmw') 
        CMD{end+1} = gMeta(preFix,file,cmd,'algo',[algo ':' ver],streamType);
    end
    for e = 1:numel(key)
        if ~isempty(strfind(preFix,':l:'))
            key = [algo ':' ver ':stage:' key{e}];
            if isempty(strfind(value{e},':$'))
                value{e} = [value{e} ':' num2str(posixtime(datetime('now','TimeZone','America/Chicago')))];
            end
        elseif ~isempty(strfind(preFix,':p:'))
            key = [algo ':' ver ':' key{e}];
        end
        CMD{end+1} = gMeta(preFix,file,cmd,key,value{e},streamType);
    end
end


