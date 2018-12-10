function [] = refreshProject(projectName)
    switch projectName
        case 'phytoSS'
            % bring newest jar from project to myJARS@phytoSS/libs
            targetPath = '/mnt/scratch1/phytoLib/myJARs/';
            refreshJava(targetPath);
            
            % bring in all jars from myJARS
            targetPath = '/mnt/scratch1/phytoLib/pAuto/algorithm/phytoSS/libs/myJARs/';    
            sourcePath = '/mnt/scratch1/phytoLib/myJARs/'
            refreshAllJava(sourcePath,targetPath);
            % add to path
            jPath = '/mnt/scratch1/phytoLib/pAuto/algorithm/phytoSS/libs/myJARs/';
            addJavaPath(jPath);
    end
end
%{
    refreshProject('phytoSS');
%}