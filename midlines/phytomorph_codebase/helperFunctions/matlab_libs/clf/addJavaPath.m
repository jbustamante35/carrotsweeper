function [] = addJavaPath(jPath)        
        %%%%%%%%%%%%%%%%%%
        % exclude list
        %%%%%%%%%%%%%%%%%%
    excludeList = {'javabuilder.jar','sqlite4java.jar','phytoG.jar'};%,'phyG-1.0-SNAPSHOT.jar'};
        %%%%%%%%%%%%%%%%%%
        % look for jars to add
        %%%%%%%%%%%%%%%%%%
    jList = gdig(jPath,{},{'jar'},1);
    %%%%%%%%%
    % add jars
    for e = 1:numel(jList)
        fprintf(['Adding@' strrep(jList{e},filesep,[filesep filesep]) '\n']);
        %%%%%%%%%%%%%%%%%%
        % parse exclude list
        %%%%%%%%%%%%%%%%%%
        toAdd = 1;        
        [p,nm,ext] = fileparts(jList{e});
        tmpName = [nm ext];
        for ex = 1:numel(excludeList)
            if strcmp(tmpName,excludeList{ex})
                toAdd = 0;
            end
        end
        %%%%%%%%%%%%%%%%%%
        % add jar
        %%%%%%%%%%%%%%%%%%
        if toAdd
            javaaddpath(jList{e});
        end
        fprintf(['Added@' strrep(jList{e},filesep,[filesep filesep]) '\n']);
    end
end

%{
jPath = '/mnt/scratch1/phytoM/phytoSS/matlab/libs/myJARs/';
addJavaPath(jPath);
%}