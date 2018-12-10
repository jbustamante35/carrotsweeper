
FilePath = '/mnt/snapper/nate/myDev/maizeWhole_mini7/output/';
FileList = {};
FileExt = {'mat'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
close all
for e = 1:100:1000
    filePath = FileList{e};
    [imagePath] = fixGraviFileName(filePath);
    imageFile = [imagePath filesep '100000.TIF'];
    imageFile2 = [imagePath filesep '000000.TIF'];
    if exist(imageFile)
        I = imread(imageFile);
        imshow(I,[])
        drawnow
    elseif exist(imageFile2)
        I = imread(imageFile2);
        imshow(I,[])
        drawnow
    else
        imageFile
    end
    
end
