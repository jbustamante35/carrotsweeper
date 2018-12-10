function [] = refreshAllJava(sourcePath,targetPath)
    %%%%%%%%%%%%%%%%    
    % remove targetPath
    rmdir(targetPath,'s');
    mkdir(targetPath);

    %%%%%%%%%%%%%%%%
    % source path
    jList = gdig(sourcePath,{},{'jar'},1);
    for e = 1:numel(jList)
        [p,nm,ext] = fileparts(jList{e});        
        targetFile = [targetPath nm];
        cmd = ['cp ' jList{e} ' ' targetFile ext];
        system(cmd);
    end
end


%{
    % refresh java is now a hetergenous call over projects
    % need to refresh over each project
    targetPath = '/mnt/scratch1/phytoLib/pAuto/algorithm/phytoSS/libs/myJARs/';    
    sourcePath = '/mnt/scratch1/phytoLib/myJARs/'
    refreshAllJava(sourcePath,targetPath);
%}
