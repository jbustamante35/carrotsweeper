FilePath = '/mnt/spaldingimages/nate/whole_TaN/';
FileList = {};
FileExt = {'tif','TIF'};
FileList = sdig(FilePath,FileList,FileExt,1);
%%
numel(FileList)
%%
I = imread(FileList{1}{1});
I = zeros([size(I) numel(FileList)]);
rm = [];
for e = 1:numel(FileList)
    try
        I(:,:,e) = imread(FileList{e}{1});
    catch
        rm = [rm e];
    end
    e
end
%%
I(:,:,rm) = [];
%%
close all
dB = cell(1,numel(FileList));
parfor e = 1:numel(FileList)
    
    tic
    tmpI = double(imread(FileList{e}{1}));
    tmpI = handleFLIP(tmpI,[]);
    tmpfI = imfilter(tmpI/255,fspecial('disk',11),'replicate');
    M = tmpfI < graythresh(tmpfI(:));
    M = bwlarge(M);
    dB{e} = bwboundaries(M);
    
    %{
    imshow(tmpI,[]);
    hold on
    plot(dB{e}{1}(:,2),dB{e}{1}(:,1),'r')
    hold off
    drawnow
    %}
    e
    toc
end
%%
cVEC = zeros(1000,2,numel(dB));
rm2 = [];
for e = 1:numel(dB)
    try
        tmpC = double(dB{e}{1});
        tmpC(tmpC(:,2) == 1,:) = [];
        sz = size(tmpC,1);
        U =  mean(tmpC([1 sz(1)],:));
        tmpC = bsxfun(@minus,tmpC,U);
        cVEC(:,:,e) = interp1(1:size(tmpC,1),tmpC,linspace(1,size(tmpC,1),1000));
    catch
        rm2 = [rm2 e];
    end
    e
end
cVEC(:,:,rm2) = [];
%%
sz = size(cVEC);
M = reshape(cVEC,[prod(sz(1:2)) sz(3)]);
[S C U E L ERR LAM] = PCA_FIT_FULL_T(M,10);
%%
close all
S = reshape(S,sz);
for e = 1:30
    plot(cVEC(:,2,e),cVEC(:,1,e),'r')
    hold on
    plot(S(:,2,e),S(:,1,e),'k')
    hold off
    axis([0 500 -100 100])
    drawnow
end
%%
 M = PCA_BKPROJ_T(C,E,U)