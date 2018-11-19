%% local focus
FilePath = '/mnt/tetra/JulianBustamante/ScottBrainard/phenotyping/NEFs/NEFs/';
FilePath = '/mnt/tetra/JulianBustamante/ScottBrainard/phenotyping/quickProcess/';
FileList = {};
FileExt = {'NEF','JPG'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% local focus
FilePath = '/mnt/tetra/JulianBustamante/ScottBrainard/phenotyping/quickProcess/NEFs/';
nFileList = {};
FileExt = {'NEF'};
nFileList = gdig(FilePath,nFileList,FileExt,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% view some images
for e = 1:10:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load some data for gmm
PER = .1;
toSample = 20000;
TS = {};
parfor e = 1:numel(FileList)
    I = imread(FileList{e});
    I = imresize(I,.1);
    I = permute(I,[3 1 2]);
    sz = size(I);
    I = reshape(I,[sz(1) prod(sz(2:3))])';
    I = I(randperm(size(I,1)),:);
    TS{e} = I(1:toSample,:);
end
%% read for dither
TILE = [];
for e = 1:10:numel(FileList)
    I = imread(FileList{e});
    TILE = [TILE ; imresize(I,.07)];
    e
end
%%
[~,map] = rgb2ind(TILE,2,'nodither');
%%
I = imread(nFileList{1});
%%
close all
Ip = rgb2ind(I,map,'nodither');
imshow(Ip,[]);
%%
MS = zeros(toSample*numel(FileList),3);
str = 1;
for e = 1:numel(TS)
    stp = str + toSample -1;
    MS(str:stp,:) = TS{e};
    str = stp + 1;
end
%% construct GMM
NC = 7;
options = statset('Display','iter','MaxIter',400);
gmm = fitgmdist(MS,NC,'Options',options);
%% apply this
I = imread(FileList{1});
Io = I;
osz = size(I);
I = permute(I,[3 1 2]);
sz = size(I);
I = reshape(I,[sz(1) prod(sz(2:3))])';
I = double(I);
kidx = gmm.cluster(I);
kidx = reshape(kidx,osz(1:2));
lRGB = label2rgb(kidx);
imshow(lRGB,[]);
%% look at overlay of clusters
NC = 7
for e = 1:NC
    out = flattenMaskOverlay(Io,kidx==e);
    imshow(out,[]);
    title(num2str(e))
    drawnow
    waitforbuttonpress
end
%%
for k = 1:3
    HI{k} = imhist(MS(:,k)/255);
end
%% get cells and sample root colors
oPath = '/mnt/tetra/JulianBustamante/ScottBrainard/phenotyping/quickReturn/';
cellClusterNumber = [0];
qrClusterLabel = 4;
blueClusterNumber = 5;
for e = 1:numel(FileList)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image and get file name information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I = imread(FileList{e});
        [pth,nm,ext] = fileparts(FileList{e});
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pre-process the image - if needed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:3
            tmp = double(I(:,:,k))/255;
            %tmp = adapthisteq(tmp);
            tmp = histeq(tmp,HI{k});
            I(:,:,k) = tmp*255;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % flip the image - if needed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if size(I,1) > size(I,2)
            I = permute(I,[2 1 3]);
            I = flip(I,2);
        end
        
        %{
        [H,T,R] = hough(cellMask,'Theta',linspace(-10,10,100));
        P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
        lines = houghlines(cellMask,T,R,P,'FillGap',100,'MinLength',100);
        I = imrotate(Io,mean([lines.theta]));
        %}

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get cell mask and cell bounding boxes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cellMask,cellBoxes] = getCellMask(I,gmm,map,cellClusterNumber);
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % process boxes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cellCrop,qrData,QRboundingBox,cropLine] = getLeftRightBox(I,cellBoxes,map,gmm,blueClusterNumber,qrClusterLabel);
        



        
        






        for r = 1:numel(cellCrop)
           
            
            cellCrop = imcrop(Io,R(r).BoundingBox);
            sub_kidx = imcrop(kidx,R(r).BoundingBox);
            [qrData,boundingBox] = getCarrotQR(cellCrop,sub_kidx,4);
            [cropLine] = getCarrotSplitLine(sub_kidx,blueClusterNumber,qrClusterLabel);


            if cropLine > boundingBox(1)
                cellCrop = flip(cellCrop,2);
                sub_kidx = flip(sub_kidx,2);
                cropLine = getCarrotSplitLine(sub_kidx,5);
                [qrData,boundingBox] = getCarrotQR(cellCrop,sub_kidx,4);
            end

            if ~isempty(qrData)



                % round one
                carrotImage = cellCrop(:,1:cropLine,:);
                carrotImage_kidx = sub_kidx(:,1:cropLine,:);
                HSV = rgb2hsv(carrotImage);
                carrotMask = HSV(:,:,2) > .2 & carrotImage_kidx ~= 2;
                carrotMask = bwlarge(carrotMask);
                iCM = sum(carrotMask,1);
                idx = find(iCM > 20);
                if ~isempty(idx)
                    cropLine = idx(end);
                    % round two to snap to edge
                    carrotImage = cellCrop(:,1:cropLine,:);
                    carrotImage_kidx = sub_kidx(:,1:cropLine,:);
                    HSV = rgb2hsv(carrotImage);
                    carrotMask = HSV(:,:,2) > .2 & carrotImage_kidx ~= 2;
                    carrotMask = bwlarge(carrotMask);
                end
                midx = find(carrotMask(:,end));
                carrotMask(midx(1):midx(end),end) = 1;
                carrotMask = imfill(carrotMask,'holes');



                try
                    [midline,contour] = getMidlineAndContour(carrotMask);
                    [vec,vecM] = sampleStraighten(midline,carrotMask,carrotImage);
                catch ME
                    fprintf(['Failed at getting, straighting, and sampling midline.\n']);
                end
                



                imshow(cellCrop,[]);
                hold on
                rectangle('Position',boundingBox,'EdgeColor','r','LineWidth',5)
                plot(cropLine*ones(size(cellCrop,1),1),1:size(cellCrop,1),'b')
                if ~isempty(midline)
                    plot(midline(:,1),midline(:,2),'g');
                    plot(contour(:,1),contour(:,2),'m');
                end
                title(qrData)
                drawnow
                hold off
                saveas(gca,[oPath qrData '.tif']);
                imwrite(vec,[oPath qrData '{imageName_straightColor}.tif']);
                imwrite(vecM,[oPath qrData '{imageName_straightMask}.tif']);
            else
                fprintf('Empty Cell.\n']);
                imshow(cellCrop,[]);
                title('Empty');
            end
            
            

        end
    catch 
    end
end
%%

%% construct GMM
NC = 3;
options = statset('Display','iter','MaxIter',400);
gmmCarrot = fitgmdist(rootStack,NC,'Options',options);


