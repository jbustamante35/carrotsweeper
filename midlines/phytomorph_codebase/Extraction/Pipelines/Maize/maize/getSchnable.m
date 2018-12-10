FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction/loganSPOOL/';
FileList = {};
FileExt = {'csv'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% 
outPath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction/schnable/';
mkdir(outPath);
for e = 1:numel(FileList)
    fidx = strfind(FileList{e},'IA');
    if ~isempty(fidx)
        fidx1 = strfind(FileList{e},'IA-');
        fidx2 = strfind(FileList{e},'_');
        iN = FileList{e}(fidx1+3:fidx2(3)-1);
        D = csvread(FileList{e});
        ridx = [find(any(abs(diff(D,1,2)) > 10*pi/180,2)) ...
                find(abs(D(:,1)) > 35*pi/180)];
        D(ridx,:) = [];
        
        csvwrite([outPath iN '.csv'],D');
        plot(D')
        hold on;
        drawnow
    end
end
%% test 
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/quickExtraction/schnable/';
FileList = {};
FileExt = {'csv'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
D = csvread(FileList{1});