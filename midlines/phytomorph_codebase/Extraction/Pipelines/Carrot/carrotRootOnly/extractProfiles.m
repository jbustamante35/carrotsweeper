function [] = extractProfiles(fileName,oPath,disp)
    [p,n,ext] = fileparts(fileName);
    I = imread(fileName);
    [QR,root] = splitBoxes(I, oPath, disp);
        
    % QRcode = readQR(QR);
    profile = phenotypeRoot(root, oPath, disp);
    % merge(QRcode, phenotype);

    
    if ~isempty(oPath)
        tmpName = [oPath n '.tif'];
        imwrite(profile, tmpName);
    end

end


%{



FilePath = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/image_preprocessing/2017/matlab_based/aligned/';
FileList = {};
FileExt = {'tif'};
FileList = gdig(FilePath,FileList,FileExt,1);
oPath = '/Users/boat/Desktop/output/profiles/';
mkdir(oPath);
toRun = numel(FileList);
for e = 1:toRun
    extractProfiles(FileList{e},oPath,false);
end


%}
        