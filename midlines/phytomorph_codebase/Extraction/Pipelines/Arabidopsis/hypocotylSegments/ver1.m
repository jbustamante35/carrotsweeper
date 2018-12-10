% scan the location for new image data
FilePath = '/mnt/spaldingdata/nate/mirror_images/forBill/rawData/';
FilePath = '/home/nate/Downloads/MOCK/';
FilePath = '/mnt/tetra/Bill_Gray/';
FilePath = [FilePath];
FileList = {};
FileExt = {'tif','TIF'};
FileList = sdig(FilePath,FileList,FileExt,1);
%FileList(1:2) = [];
%%
close all
oPath = '/mnt/spaldingdata/nate/mirror_images/forBill/return/';
oPath = '~/Downloads/returnForBill/';
oPath = '/mnt/tetra/forBillReturn/';
FilePath = '/home/nate/Downloads/returnMOCK/';
FilePath = '/mnt/tetra/Bill_Gray/';
mkdir(oPath);
for f = 1:numel(FileList)
    op0(FileList{f},oPath,1);
    f
end
