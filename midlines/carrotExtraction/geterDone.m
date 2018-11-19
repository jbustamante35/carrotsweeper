%% local focus
FilePath = '/mnt/tetra/JulianBustamante/carrots/phenotyping/quickProcess/NEFs/NEFs/';
FileList = {};
FileExt  = {'NEF'};
FileList = gdig(FilePath,FileList,FileExt,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% view some images
for e = 1:10:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load some data for gmm
PER = .1;
toSample = 20000;
TS = {};
parfor e = 1:numel(FileList)
    I  = imread(FileList{e});
    I  = imresize(I,.1);
    I  = permute(I,[3 1 2]);
    sz = size(I);
    I  = reshape(I,[sz(1) prod(sz(2:3))])';
    I  = I(randperm(size(I,1)),:);
    TS{e} = I(1:toSample,:);
end

%% 
MS = zeros(toSample*numel(FileList),3);
str = 1;
for e = 1:numel(TS)
    stp = str + toSample -1;
    MS(str:stp,:) = TS{e};
    str = stp + 1;
end

%% construct GMM
NC      = 7;
options = statset('Display','iter','MaxIter',400);
gmm     = fitgmdist(MS,NC,'Options',options);

%% apply this
I    = imread(FileList{1});
Io   = I;
osz  = size(I);
I    = permute(I,[3 1 2]);
sz   = size(I);
I    = reshape(I,[sz(1) prod(sz(2:3))])';
I    = double(I);
kidx = gmm.cluster(I);
kidx = reshape(kidx,osz(1:2));
lRGB = label2rgb(kidx);
imshow(lRGB,[]);

%% look at overlay of clusters
NC = 7;
for e = 1:NC
    out = flattenMaskOverlay(Io,kidx==e);
    imshow(out,[]);
    title(num2str(e));
    drawnow;
    waitforbuttonpress;
end

%% get cells and sample root colors
for e = e:numel(FileList)
    I = imread(FileList{e});
    
    
    

    Io = I;
    osz = size(I);
    I = permute(I,[3 1 2]);
    sz = size(I);
    I = reshape(I,[sz(1) prod(sz(2:3))])';
    I = double(I);
    kidx = gmm.cluster(I);
    kidx = reshape(kidx,osz(1:2));
    
    
    
    cellClusterNumber = 2;
    borderClusterNumber = 2;
    smallThreshold = 500000;
    cellMask = kidx == cellClusterNumber;
    cellMask = imclose(cellMask,strel('square',51));
    cellMask = bwareaopen(cellMask,smallThreshold);
    
    [H,T,R] = hough(cellMask,'Theta',linspace(-10,10,100));
    P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
    lines = houghlines(cellMask,T,R,P,'FillGap',100,'MinLength',100);
    I = imrotate(Io,mean([lines.theta]));
    
    
    
    
    
    Io = I;
    osz = size(I);
    I = permute(I,[3 1 2]);
    sz = size(I);
    I = reshape(I,[sz(1) prod(sz(2:3))])';
    I = double(I);
    kidx = gmm.cluster(I);
    kidx = reshape(kidx,osz(1:2));
    
    
    
    
    
    
    imshow(cellMask,[]);
    
    
    
    R = regionprops(logical(cellMask),'BoundingBox');
    rootStack = [];
    
    
    for r = 1:numel(R)
        %{
        imshow(Io,[])
        hold on
        rectangle('Position',R(r).BoundingBox);
        hold off
        %}
        cellCrop = imcrop(Io,R(r).BoundingBox);
        sub_kidx = imcrop(kidx,R(r).BoundingBox);
        [qrData,boundingBox] = getCarrotQR(cellCrop,sub_kidx,4);
        [cropLine] = getCarrotSplitLine(sub_kidx,5);
        if ~isempty(qrData)



            % round one
            carrotImage = cellCrop(:,1:cropLine,:);
            carrotImage_kidx = sub_kidx(:,1:cropLine,:);
            HSV = rgb2hsv(carrotImage);
            carrotMask = HSV(:,:,2) > .2 & carrotImage_kidx ~= 2;
            carrotMask = bwlarge(carrotMask);
            iCM = sum(carrotMask,1);
            idx = find(iCM > 60);
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

            %{
            csz = size(carrotImage);
            carrotImageP = permute(carrotImage,[3 1 2]);
            sz = size(carrotImageP);
            carrotImageP = reshape(carrotImageP,[sz(1) prod(sz(2:3))])';
            carrotImageP = double(carrotImageP);
            rootStack = [rootStack;carrotImageP];
            %}


            try
                midline = [];
                contour = [];
                [out] = isolate_carrot_Roots(double(~carrotMask),0,[],[]);
                midline = out(1).midlines.data';
                sz = size(carrotMask);
                midline(:,1) = midline(:,1) - sz(2)/2;
                midline(:,1) = -midline(:,1);
                midline(:,1) = midline(:,1) + sz(2)/2;

                contour = out(1).contours.data';
                contour(:,1) = contour(:,1) - sz(2)/2;
                contour(:,1) = -contour(:,1);
                contour(:,1) = contour(:,1) + sz(2)/2;


                [DomainS,DomainG] = extendCarrotMidline(midline,[0 0],carrotMask);
                dsz = size(DomainG);
                vec = [];
                for k = 1:size(carrotImage,3)
                    vec(:,k) = ba_interp2(double(carrotImage(:,:,k))/255,DomainS(:,2),DomainS(:,1));
                end

                for k = 1:size(carrotImage,3)
                    vec(:,k) = ba_interp2(double(carrotImage(:,:,k))/255,DomainS(:,2),DomainS(:,1));
                end


                vec = reshape(vec,[dsz(1) dsz(2) 3]);
            catch ME
                ME;
            end
            %waitforbuttonpress

            %{

            if exist('gmmCarrot')
                cidx = gmmCarrot.cluster(carrotImageP);
                cidx = reshape(cidx,csz(1:2));
                out = flattenMaskOverlay(carrotImage,cidx==0);
            end
            %}
            
            
            
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
            
        end
            
    end
end

%%

%% construct GMM
NC = 3;
options = statset('Display','iter','MaxIter',400);
gmmCarrot = fitgmdist(rootStack,NC,'Options',options);


