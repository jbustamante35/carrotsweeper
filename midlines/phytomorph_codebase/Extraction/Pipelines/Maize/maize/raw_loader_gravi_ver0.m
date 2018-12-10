function [] = raw_loader_gravi_ver0()
    %% find the mat files to load from
    FilePath = '/mnt/snapper/nate/myDev/maizeWhole_mini7/';
    FileList = {};
    FileExt = {'mat'};
    verbose = 1;
    FileList = gdig(FilePath,FileList,FileExt,verbose);
    %%%%%%%%%%%%%%%%%%%%%
    % for each file
    %%%%%%%%%%%%%%%%%%%%%
    parfor m = 1:numel(FileList)
        fprintf(['Starting on load of:' num2str(m) '\n']);
        ERR(m) = graviExtractCore(FileList{m},m);
    end
end