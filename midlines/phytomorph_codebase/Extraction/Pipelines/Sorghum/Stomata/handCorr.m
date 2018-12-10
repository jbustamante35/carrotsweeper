%%
massDownload('/iplant/home/leakey_cyverse/quickReturn_sorghum2016_verFinal11/', 'locations.csv','/home/nate/Download/testB/');
%%
Data = readtext('/home/nate/Downloads/CountTest_sorghum.csv')
dPath = '/home/nate/Download/testB/';
cnt = 1;
HAND = [];
AUTO = [];
for e = 2:size(Data,1)
    f1 = Data{e,1};
    f2 = Data{e,2};
    f3 = Data{e,3};
    for r = 1:(4-length(f2))
        f2 = ['0' f2];
    end
    fileNameT = Data{e,5};
    f2 = [fileNameT '_WITH2locations.csv'];
    f1 = [fileNameT '_locations.csv'];
    try
    tmp2 = csvread([dPath f2]);
    tmp1 = csvread([dPath f1]);
    HAND(cnt) = Data{e,6};
    AUTO(cnt,1) = size(tmp1,1);
    AUTO(cnt,2) = size(tmp2,1);
    KL{cnt} = f2;
    cnt = cnt+1;
    catch
    end
end
corr(HAND',AUTO)
%%
[J,idx] = sort(abs(AUTO(:,1) - HAND'),'descend');
I = imread(['/iplant/home/leakey_cyverse/quickReturn_sorghum2016_verFinal10/' strrep(KL{idx(3)},'WITH2locations.csv','labeledImage.jpg')]);
close all
imshow(I,[]);
%%
close all
plot(HAND,AUTO(:,1),'.');
%%
FilePath = '/home/nate/Download/testIIII/';
cFileList = {};
FileExt = {'csv'};
cFileList = gdig(FilePath,cFileList,FileExt,1);
%%
for e = 1:numel(cFileList)
    [cp{e} cnm{e}] = fileparts(cFileList{e});
    fidx = strfind(cnm{e},'_');
    cnm{e} = cnm{e}(1:fidx(end)-1);
end