function [] = localClean(pthList)
    %%%%%%%%%%%%%%%%%%
    % clean local
    for ele = 1:numel(pthList)
        %%%%%%%%%%%%%%%%%%
        % store old path    
        p = pwd;
        %%%%%%%%%%%%%%%%%%
        rmdir(pthList{ele},'s');
        
        %if ~(isempty(find(pthList{ele},'working')))
        %    mkdir(pthList{ele})
        %end
        %{
        % rm directories and files - force
        cmd = ['rm -f -r ' pthList{ele} '*'];        
        %fprintf(['eval@' cmd '\n']);
        system(cmd);        
        % restore path
        cmd = ['cd ' p];
        [status,resul] = system(cmd);
        %}
    end
end