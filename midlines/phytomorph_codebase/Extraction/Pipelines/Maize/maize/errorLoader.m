function [] = errorLoader()
    % define errorpath
    ePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/massExtraction/errorPath/';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% scan for new images
    FileList = {};
    FileExt = {'csv'};
    verbose = 1;
    FileList = gdig(ePath,FileList,FileExt,verbose);
    for e = 1:numel(FileList)
        fidx = strfind(FileList{e},'~');
        pathString = strrep(FileList{e}(fidx(end)+1:end),'--','/');
        pathString(end-3:end) = [];
        pathString = [pathString '/'];
        main(pathString);
        errorString = FileList{e}(fidx(1)+1:fidx(2)-1)
    end
end