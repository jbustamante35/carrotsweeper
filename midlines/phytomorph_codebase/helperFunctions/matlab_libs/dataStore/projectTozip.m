function [] = projectTozip(fileList,zipFile)
    % import libs
    import phytoG.locked.Bpersist.Bfs.implementations.*;
    % create file list
    fileList = toBfile(fileList);
    % create file system
    localFS = Bfs_local();
    % zip
    localFS.zipp(fileList, zipFile);
end
