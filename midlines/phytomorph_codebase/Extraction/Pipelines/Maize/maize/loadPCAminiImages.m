% define inpath 
inFilePath = '/mnt/spaldingdata/nate/mirror_images/beta_miniMaize/';
%%% scan for new images
FileList = {};
FileExt = {'tiff','TIF','tif'};
verbose = 1;
SET = gdig(inFilePath,FileList,FileExt,verbose);
I = imread(SET{1});
I = imresize(I,.5);
[I BOX] = imcrop(I);
STACK = zeros([size(I) numel(SET)]);
rm = logical(zeros(1,numel(SET)));
for e = 1:numel(SET)
    try
        tmp = imread(SET{e});
        tmp = imresize(tmp,.5);
        tmp = imcrop(tmp,BOX);
        STACK(:,:,e) = tmp;
        e
        numel(SET)
    catch
        rm(e) = 1;
    end
end
STACK(:,:,rm) = [];
%%
STACK = reshape(STACK,[size(STACK,1)*size(STACK,2) size(STACK,3)]);
[S C U E L ERR LAM] = PCA_FIT_FULL(STACK',200);
S = reshape(S',[size(tmp,1) size(tmp,2) size(S,1)]);