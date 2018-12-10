function [] = refreshJava(targetPath)
    %%%%%%%%%%%%%%%%
    % source path
    sourceFile = '/mnt/scratch1/phytoLib/myJavaBuilds/phytoGraph/dist/phytoG.jar';
    [pth,nm,ext] = fileparts(sourceFile);    

    %%%%%%%%%%%%%%%%
    % remove and restore    
    targetFile = [targetPath nm ext];
    mkdir(targetPath);
    
    %%%%%%%%%%%%%%%%
    % remove and restore
    delete(targetFile);
    cmd = ['cp ' sourceFile ' ' targetFile];
    system(cmd);
    
end

%{
    % refresh java is now a hetergenous call over projects
    % need to refresh over each project
    targetPath = '/mnt/spaldingdata/nate/phytoLib/pAuto/algorithm/phytoSS/libs/myJARs/';
    refreshJava(targetPath);
%}
