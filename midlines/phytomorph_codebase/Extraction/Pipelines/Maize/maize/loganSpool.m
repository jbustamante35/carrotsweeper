%% dig
FilePath = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/allExtraction/csvSPOOL/';
FileList = {};
FileExt = {'csv'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% open and resave for logan
loganOut = '/mnt/spaldingdata/Takeshi/allMaizeMovies_results/allExtraction/loganSPOOL/';
mkdir(loganOut);
close all
for i = 1:numel(FileList)
    try
        [p n ext] = fileparts(FileList{i});
        % snip to folder
        fidx = strfind(n,'----');
        newName = n(fidx(end) + 4:end);
        fidx = strfind(newName,'_');
        plateName = newName(1:fidx(2)-1);
        wellNames = newName(fidx(2)+1:end);
        D = csvread(FileList{i});
        newFileName = [loganOut newName(1:fidx(end)-1) '.txt'];
        csvwrite(newFileName,D);
        plot(D');
        drawnow
    catch
        
    end
end
