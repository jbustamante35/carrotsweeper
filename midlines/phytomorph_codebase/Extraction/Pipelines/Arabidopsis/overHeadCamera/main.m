%% whole run
FilePath = '/mnt/tetra/nate/overHead/';
sFileList = {};
FileExt = {'jpg'};
sFileList = sdig(FilePath,sFileList,FileExt,1);
%% show movie
s = 30;
cnt = 1;
M = [];
for e = 1:10:numel(sFileList{s})
    I = imread(sFileList{s}{e});
    if ~isNight(I,50)
        M(:,:,:,cnt) = I;
        cnt = cnt + 1;
        imshow(I,[]);
    end
end
%%
close all
for e = 1:size(M,4)
    imshow(M(:,:,:,e)/255,[]);
    drawnow
end
%%
FilePath = '/mnt/tetra/nate/overHead/';
FileList = {};
FileExt = {'jpg'};
FileList = gdig(FilePath,FileList,FileExt,1);
%% sample random images
close all
mag = .2;
FileList = FileList(randperm(numel(FileList)));
S = [];
cnt = 1;
e = 1;
TILE = [];
while cnt < 100
    I = imread(FileList{e});
    I = imresize(I,mag);
    aI = mean(I(:));
    e = e + 1;
    if aI > 40
        sz = size(I);
        imshow(I,[]);
        TILE = [TILE;I];
        drawnow
        I = reshape(I,[prod(sz(1:2)) sz(3)]);
        S = [S;I];
        cnt = cnt + 1;
        cnt
    end
end
S = double(S);
%% grab checkboards
close all
cS = [];
checkBoardList = {};
for e = 1:numel(sFileList)
    I = imread(sFileList{e}{1});
    checkBoardList{end+1} = sFileList{e}{2};
    I = imresize(I,mag);
    aI = mean(I(:));
    if aI > 40
        imshow(I,[]);
        TILE = [TILE;I];
        title(sFileList{e}{1});
        drawnow
        sz = size(I);
        I = reshape(I,[prod(sz(1:2)) sz(3)]);
        cS = [cS;I];
        cnt = cnt + 1;
    end
end
cS = double(cS);
%% stack images from random and checker boards
mS = [S;cS];
%% geneate map for colors
[X,map] = rgb2ind(double(TILE)/255,9);
%% try main
oPath = '/mnt/tetra/nate/forCullenMovie/';

FilePath = '/mnt/tetra/nate/AlgoTest/';
FilePath = '/mnt/tetra/nate/arco/';
sFileList = {};
FileExt = {'jpg'};
sFileList = sdig(FilePath,sFileList,FileExt,1);
%%
overHead_main(sFileList{1},oPath,map,GMModel);
%% publish
overHeadFunc = @(X)overHead_main(X,'./output/',map,GMModel);
pF = partialFunction(overHeadFunc,'overHeadFuncAPP2');
pF.publish();
%% try mp4
V = VideoReader('/mnt/tetra/nate/timelapse-20181029-09-36-08.mp4');
%% find the data for labels and crop boxes
parfor stack = 51:60
    [typeTable{stack},cropTableCheckerBoard{stack},cropTableRedSquares{stack}] = parseAndLabelStack(sFileList{stack},map,false);
end
%%
for stack = 51:60%1:numel(typeTable)
    try
        % get the file names for the day images only
        fileNames = getDayfileNames(typeTable{stack});
        [pth,nm,ext] = fileparts(fileNames{1});
        fidx = strfind(nm,'_');
        nm = nm(1:fidx(1)-1);
        fidx = strfind(pth,filesep);
        pth = pth(fidx(end-1)+1:fidx(end)-1);
        % get the names of the checker board images
        CheckerBoardfileNames = getCheckBoardfileNames(typeTable{stack});
        % get the crop boxes for the checker board images
        checkerBoardBoxes = getFieldForFileName(CheckerBoardfileNames,cropTableCheckerBoard{stack},'CropBoxes');
        % 
        [tform,resEstimate] = checkerBoardAnalysis(CheckerBoardfileNames{1},checkerBoardBoxes{1},true);
        % get crop boxes for the first N image
        N = 8;
        cropBOX = {};
        rm = [];
        for f = 1:N
            [cropBOX{f}] = getPlantBoxesForFile(cropTableRedSquares{stack},fileNames{f});
            if numel(cropBOX{f}) ~= 9
                rm(f) = true;
            else
                rm(f) = false;
            end
        end
        cropBOX(find(rm)) = [];
        [cropBOX] = mergeCropBoxes(cropBOX);
        renderBoxesOnImage(fileNames{6},cropBOX,'r',3);
        outTable = table;
        Area = [];
        for f = 1:numel(fileNames)
            I = double(imread(fileNames{f}))/255;
            [Area(f,:),QR_text{f}] = getPlantCell_and_Mask(I,cropBOX,GMModel,map);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather possible fields
        FLDS = {};
        for f = 1:numel(QR_text)
            for b = 1:numel(QR_text{f})
                [flds] = mj_getFields(QR_text{f}{b});
                FLDS = [FLDS flds];
                FLDS = unique(FLDS);
            end
        end
        % get values
        HEADER = FLDS';
        for f = 1:numel(QR_text)
            for b = 1:numel(QR_text{f})
                for h = 1:size(HEADER,1)
                    [~,value] = mj_getFields(QR_text{f}{b},HEADER{h,1});
                    if b == 1 | ~isempty(value{1})
                        HEADER{h,b+1} = value{1};
                    end
                end                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        


        oArea = {};
        for tm = 1:size(Area,1)
            for p = 1:numel(cropBOX)
                oArea{p}(tm) = Area(tm,p);
            end
        end
        outTable = table(fileNames,oArea{1}(:),oArea{2}(:),oArea{3}(:),oArea{4}(:),oArea{5}(:),oArea{6}(:),oArea{7}(:),oArea{8}(:),oArea{9}(:));
        outTable = table2cell(outTable);
    
        writetable(outTable,[oPath pth '__' nm '.csv']);
    catch ME
        stack
        ME
    end
end

%% crop out the red squares, remove white and sample colors
toCollect = false;
if toCollect
    sCS = [];
end
for e = 1:numel(cropTableRedSquares)
    rnd = randi(size(cropTableRedSquares{e}.fileName,1),10,1);
    for r = 1:10
        fn = cropTableRedSquares{e}.fileName{rnd(r)};
        cropBOX = cropTableRedSquares{e}.CropBoxes{rnd(r)};
        cropBOX(1:2) = cropBOX(1:2) + 20;
        cropBOX(3:4) = cropBOX(3:4) - 40;
        I = double(imread(fn))/255;
        subI = imcrop(I,cropBOX);
        if ~toCollect
            
            sz = size(subI);
            cidx = GMModel.cluster(reshape(subI,[prod(sz(1:2)) sz(3)]));
            subL = reshape(cidx,sz(1:2));
            
            
            plantMask = subL == 3;
            plantMask = bwareaopen(plantMask,300);
            plantMask = imclearborder(plantMask);
            %{
            midx = find(plantMask);
            sig = [];
            for k = 1:3
                tmp = subI(:,:,k);
                sig = [sig tmp(midx)];
            end
            cidx = GMModel2.cluster(sig);
            plantMask = zeros(size(plantMask));
            plantMask(midx) = cidx;
            
            out = label2rgb(plantMask);
            %}
            out = flattenMaskOverlay(subI,plantMask);
            
            
            %subL = rgb2ind(subI,newMap,'nodither');
            %out = label2rgb(subL);
            
            
            
            imshow(out,[]);
            drawnow
        end
        
        if toCollect
            subI = imresize(subI,[200 200]);
            sCS = [sCS subI];
        end
        
    end
    e
end
%%
sz = size(sCS);
sCSr = reshape(sCS,[prod(sz(1:2)) sz(3)]);
GMModel = fitgmdist(sCSr(1:10:end,:),4);
gidx = GMModel.cluster(sCSr);
GMModel2 = fitgmdist(sCSr(gidx==3,:),3);
%% make map - new
[~,newMap] = rgb2ind(sCS,7,'nodither');
%%
oPath = '/mnt/tetra/nate/forCullenMovie/';
for stack = 1:numel(sFileList)
    [pth,nm,ext] = fileparts(sFileList{stack}{1});
    fidx = strfind(nm,'_');
    nm = nm(1:fidx(1)-1);
    vFile = [oPath num2str(nm) '.avi'];
    
    v = VideoWriter(vFile);
    open(v);
    
    
    for f = 1:numel(sFileList{stack})
        I = double(imread(sFileList{stack}{f}))/255;
        if ~isNight(I,50)
            writeVideo(v,I);
        end
        f
    end

    close(v)
    stack
end


%%
close all
for s = 37:numel(sFileList)
    
    I = double(imread(sFileList{s}{10}))/255;
    for r = 1:4
        I(1,:,1) = 1;
        I(1,:,2) = 0;
        I(1,:,3) = 0;
        I = imrotate(I,90);
    end
    L_nodither = rgb2ind(I,map,'nodither');
    L_dither = rgb2ind(I,map,'dither');
    
    RGB = label2rgb(L);
    
    redTape = L_dither == 4 | L_dither == 6 | L_dither == 2;
    redTape = bwlarge(redTape);
    redTape = imclose(redTape,strel('square',71));
    
    cropBOXES = imfill(redTape,'holes') - redTape;
    cropBOXES = bwareaopen(cropBOXES,250000);
    cropBOXES = bwlarge(cropBOXES,12);
    cropBOXES_Region = regionprops(logical(cropBOXES),'BoundingBox','Area');
    
    cbArea = [];
    for i = 1:numel(cropBOXES_Region)
        cbArea(i) = cropBOXES_Region(i).BoundingBox(3)*cropBOXES_Region(i).BoundingBox(4);
        cbmaxDIM(i) = max(cropBOXES_Region(i).BoundingBox(3:4));
    end
    rm = cbArea > 700000 | cbmaxDIM > 800;
    cropBOXES_Region(rm) = [];
    
    
    
    
    blueSquare = L_dither == 7 | L_dither == 8;
    blueSquare = bwlarge(blueSquare);
    blueSquare = imclose(blueSquare,strel('square',71));
    blueSquare = imclearborder(blueSquare);
    
    
    innerBlueSquare = imfill(blueSquare,'holes') - blueSquare;
    innerBlueSquare = bwlarge(innerBlueSquare,1);
    innerBlueSquare = bwareaopen(innerBlueSquare,10000);
    innerBlueSquare_Region =  regionprops(logical(innerBlueSquare),'BoundingBox');
    
    
    
    white = L_dither == 3;
    whiteP = bwareaopen(white,1000);
    whiteP = imclose(whiteP,strel('square',51));
    whiteP = imclearborder(whiteP);
    
    other = ~white & ~blueSquare & ~redTape;
    
    %{
    out = flattenMaskOverlay(I,redTape,1,'r');
    %out = flattenMaskOverlay(out,cropBOXES,.8,'m');
    
    out = flattenMaskOverlay(out,blueSquare,1,'b');
    out = flattenMaskOverlay(out,white,1,'w');
    out = flattenMaskOverlay(out,other,.1,'m');
    %}
    
    imshow(I,[]);
    hold on
    if isempty(innerBlueSquare_Region)
        
        for e = 1:numel(cropBOXES_Region)
            rectangle('Position',cropBOXES_Region(e).BoundingBox,'EdgeColor','g','LineWidth',3);
        end
        
    else
        for e = 1:numel(innerBlueSquare_Region)
            rectangle('Position',innerBlueSquare_Region(e).BoundingBox,'EdgeColor','c','LineWidth',3);
        end
    end
    hold off
    drawnow
end