function [] = mainManyCarrot(fileName,oPath,disp)
    [p,n,ext] = fileparts(fileName);
    wholeImage = imread(fileName);
    
    % Generate cell array with cropped and aligned carrots
    
    boxes = extractBoxes(wholeImage,oPath,disp);
    
    for i = 1:numel(boxes)
        
        [QR,root] = splitBoxes(boxes{i},oPath,disp);
        
        % QRcode = readQR(QR);
        profile = phenotypeRoot(root, oPath, disp)
        % merge(QRcode, profile);
        
%         if ~isempty(oPath)
%             tmpNameCell = [oPath n '_' num2str(i) '.tif'];
%             imwrite(boxes{1}, tmpName);
%         end
    end

close all

end


%{


FilePath = '/Users/boat/Dropbox/Wisconsin/Thesis/Phenotyping/Raws/2017/NEFs/training/';
FileList = {};
FileExt = {'NEF'};
FileList = gdig(FilePath,FileList,FileExt,1);
oPath = '/Users/boat/Desktop/output/';
mkdir(oPath);
toRun = numel(FileList);
toRun = 1
for e = 1:toRun
    mainManyCarrot(FileList{e},oPath,true);
end

%}
        