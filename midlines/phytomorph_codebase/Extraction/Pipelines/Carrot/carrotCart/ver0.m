FilePath = '/home/nate/forCarrotCart/';
FileList = {};
FileExt = {'JPG'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
I = imread(FileList{1});
%%
[J,BOX] = imcrop(I);
%%
clear I
for e = 1:numel(FileList)
    tmp = imread(FileList{e});
    I(:,:,:,e) = imcrop(tmp,BOX);
    imshow(I(:,:,:,e),[])
    drawnow
end
%%
for e = 1:size(I,4)
    mI(:,:,:,e) = imresize(I(:,:,:,e),.1);
end
%%
mI = permute(mI,[1 2 4 3]);
sz = size(mI);
mI = reshape(mI,[prod(sz(1:3)) sz(4)]);
%%
options = statset('Display','Iter');
SKIP = 10;
obj = gmdistribution.fit(double(mI(1:SKIP:end,:)),2,'Options',options);
kidx = obj.cluster(double(mI));
sidx = find(kidx==1);
obj2 = gmdistribution.fit(double(mI(sidx,:)),2,'Options',options);

%%
oPath = '/home/nate/quickCart/';
mkdir(oPath);
for e = 1:size(I,4)
    tmp = I(:,:,:,e);
    sz = size(tmp);
    tmp = reshape(tmp,[prod(sz(1:2)) sz(3)]);
    kidx = obj.cluster(double(tmp));
    sidx = find(kidx==1);
    kidx2 = obj2.cluster(double(tmp(sidx,:)));
    kidx(sidx(kidx2==1)) = 0;
    kidx = reshape(kidx+1,sz(1:2));
    %kidx = label2rgb(kidx);
    %imshow(kidx,[]);
    out = flattenMaskOverlay(double(I(:,:,:,e))/255,kidx==1,.5,'r');
    imshow(out,[])
    drawnow
    imwrite(out,[oPath num2str(e) '_overlay.tif']);
    imwrite(double(I(:,:,:,e))/255,[oPath num2str(e) '.tif']);
    e
end