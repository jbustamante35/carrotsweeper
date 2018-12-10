function [CMD] = gMeta(pre,file,cmd,key,value,baseSwitch,nvalue)
    % base cmd
    switch baseSwitch
        case 'single'
            CMD = 'imeta #cmd# -d "#file#" "#key#" "#value#"';
        case 'stream'
            CMD = '#cmd# -d "#file#" "#key#" "#value#"';
            %CMD = ['#key#=#value#'];
    end
    
    % 
    %if ~isempty(strfind(value,'TM'))
        %{
        fidx = strfind(value,':');
        TMCMD = 'TM=$(date +%Y_%m_%d_%H_%M_%S)';
        CMD = strrep(CMD,'"#value#"','#value#');
        CMD = strrep(CMD,'#value#',['$TM' ':' value((fidx(1)+1):end)]);
        
        switch baseSwitch
            case 'single'
                CMD = [TMCMD ';' CMD];
            case 'stream'
                CMD = strrep(
        end
        %}
    if strcmp(value,'%')
        CMD = strrep(CMD,'"#value#"','#value#');
        CMD = strrep(CMD,'#value#',value);
    else
        CMD = strrep(CMD,'#value#',value);
    end
    
    
    
    if nargin == 7
        CMD = [CMD ' #nvalue#'];
        CMD = strrep(CMD,'#nvalue#',nvalue);
    end
    
    
    
    CMD = strrep(CMD,'#cmd#',cmd);
    CMD = strrep(CMD,'#file#',file);
    CMD = strrep(CMD,'#key#',[pre key]);
    
end