FilePath = '//mnt/scratch4/For Nate_Temp/';
FileExt = {'JPG'};
FileList = {};
verbose = 1;
[FileList] = gdig(FilePath,FileList,FileExt,verbose);
kidx = [];
for e = 1:numel(FileList)
    if ~isempty(findstr(FileList{e},'mask'))
        kidx = [kidx e];
    end
end
FileList(kidx) = [];
%%
CL = [];
for e = 1:numel(FileList)
    tmp = imread(FileList{e});
    tmp = imresize(tmp,.5);
    tmp = reshape(tmp,[size(tmp,1)*size(tmp,2) size(tmp,3)]);
    CL = [CL;tmp];
    e
end
%%
sam = 100;
options = statset('Display','iter');
gm = fitgmdist(double(CL(1:sam:end,:)),5,'Options',options);

%%
for e = 1:numel(FileList)
    I = imread(FileList{e});   
    sz = size(I);
    tmp = reshape(I,[sz(1)*sz(2) sz(3)]);
    tmp = cluster(gm,double(tmp));
    tmp = reshape(tmp,sz(1:2));
    for k = 1:5
        imshow(tmp==k,[])
        waitforbuttonpress
    end
    tmp = label2rgb(tmp);
    
end
%%
for e = 1:numel(FileList)
    I = imread(FileList{e});   
    I = imcrop(I);
    sz = size(I);
    tmp = reshape(I,[sz(1)*sz(2) sz(3)]);
    tmp = cluster(gm,double(tmp));
    tmp = reshape(tmp,sz(1:2));        
    imshow(tmp==1,[]);
    
    
    
    tmp = label2rgb(tmp);
    
    
end
%%
kidx = cluster(gm,double(CL));
fidx = find(kidx == 1);
gm2 = fitgmdist(double(CL(fidx(1:10:end),:)),2,'Options',options);
%%
for e = 1:numel(FileList)
    I = imread(FileList{e});   
    sz = size(I);
    tmp = reshape(I,[sz(1)*sz(2) sz(3)]);
    tmpK = cluster(gm,double(tmp));
    fidx = find(tmpK==1);
    tmp2 = cluster(gm2,double(tmp(fidx,:)));
    tmpK(fidx) = tmp2 + 5;
    tmpK = reshape(tmpK,sz(1:2));
    %{
    for k = 1:5
        imshow(tmp==k,[])
        waitforbuttonpress
    end
    %}
    tmpC = label2rgb(tmpK);
    imshow(tmpC,[]);
    
end





