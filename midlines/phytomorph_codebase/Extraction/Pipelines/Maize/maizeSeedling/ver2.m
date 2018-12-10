% scan for image files
FilePath = '//mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/photos/';
FilePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/photos/growthchamber12.7.15/';
FilePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/seedlingData/growthchamber12.7.15/';

FileList = {};
FileExt = {'jpeg','tif'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% view images
for e = 1:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow
end
%% loop over images
oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/return/seedlings/output/';
for e = 1:numel(FileList)
    singleSeedlingImage(FileList{e},50,5,100,100,4,oPath);
end