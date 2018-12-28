%% local focus
FilePath = '/mnt/tetra/JulianBustamante/ScottBrainard/phenotyping/NEFs/NEFs/';
FilePath = '/mnt/tetra/JulianBustamante/ScottBrainard/phenotyping/quickProcess/';
FilePath = '/mnt/tetra/JulianBustamante/carrots/phenotyping/quickProcess/';
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
oPathFrame = '/mnt/tetra/JulianBustamante/ScottBrainard/phenotyping/quickReturn2/';

oPath = '/mnt/tetra/JulianBustamante/carrots/phenotyping/quickReturn/';
oPathFrame = '/mnt/tetra/JulianBustamante/carrots/phenotyping/quickReturn2/';
mkdir(oPathFrame)
cellClusterNumber = [0];
qrClusterLabel = 4;
blueClusterNumber = 5;
for e = 1:numel(FileList)
    try
        
        filePushList = {};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % read the image and get file name information
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        I = imread(FileList{e});
        [pth,nm,ext] = fileparts(FileList{e});
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % frame builder
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [OUT_SQUARE,IN_SQUARE,RED_SQUARE,RIGHT_STRIP_SQUARE,HEADER_SQUARE,CARROT_SQUARE,ORIN] = findDarkLines(I,e,oPathFrame);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get qr msg
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [qrMSG] = getQRmsg(I,RED_SQUARE);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % straighten and sample
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [R,carrotImage,carrotMask,midline,contour,straightRGB,straightMSK] = cropTraceStraightenSample(I,CARROT_SQUARE,ORIN,qrMSG);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % new JSON string format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % strip ticket
        linkTo = stripiTicket(FileList{e});
        %linkPath = stripiTicket(rPath);
        [jP,tN] = fileparts(FileList{e});
        N = 4000;
        for b = 1:numel(qrMSG)
            
            if ~isempty(qrMSG{b})
                
                widthProfile = sum(flip(straightMSK{b},1),2);
                widthProfile = [widthProfile;zeros(N-numel(widthProfile),1)];
                
                area = sum(carrotMask{b}(:));
                
                
                % generate json style format document
                phenoTypeDocument = [];
                phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,linkTo,'orginalImage','orginalImage');
                phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,widthProfile,'rootWidthProfile','rootWidthProfile');
                phenoTypeDocument = generatePhenotypeNode(phenoTypeDocument,area,'rootArea','rootArea');

                % generate json string
                JSON_string = savejson('carrotDoc',phenoTypeDocument);

                % save json document
                filePushList{end+1} = [oPath qrMSG{b} '{output_json}.json'];
                fileID = fopen(filePushList{end},'w');
                fprintf(fileID,strrep(JSON_string,'\/','\\/'));
                fclose(fileID);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % new JSON string format
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    catch ME
        ME
    end
end
       
%%
%{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pre-process the image - if needed
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for k = 1:3
            tmp = double(I(:,:,k))/255;
            tmp = adapthisteq(tmp);
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get cell mask and cell bounding boxes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cellMask,cellBoxes] = getCellMask(I,gmm,map,cellClusterNumber);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rotate the image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [H,T,R] = hough(cellMask,'Theta',linspace(-10,10,100));
        P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
        lines = houghlines(cellMask,T,R,P,'FillGap',100,'MinLength',100);
        I = imrotate(I,mean([lines.theta]));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get cell mask and cell bounding boxes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cellMask,cellBoxes] = getCellMask(I,gmm,map,cellClusterNumber);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % process boxes
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cellCrop,sub_kidx,qrData,QRboundingBox,cropLine] = getLeftRightBox(I,cellBoxes,map,gmm,blueClusterNumber,qrClusterLabel);
        
        rootStack = [];
        carrotImage = {};
        carrotMask = {};
        for r = 1:numel(cellCrop)
            if ~isempty(qrData{r})
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get carrot mask and carrot image
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [carrotMask{r},carrotImage{r}] = getCarrotRootMask(cellCrop{r},sub_kidx{r},cropLine(r));

                try
                    midline = [];
                    contour = [];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % get midline and contour
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [midline,contour] = getMidlineAndContour(carrotMask{r});
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % sample straight midline
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [vec,vecM] = sampleStraighten(midline,carrotMask{r},carrotImage{r});
                catch ME
                     fprintf(['Failed at getting, straighting, and sampling midline.\n']);
                end


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % start display
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                imshow(cellCrop{r}/255,[]);
                hold on
                rectangle('Position',QRboundingBox{r},'EdgeColor','r','LineWidth',5)
                plot(cropLine(r)*ones(size(cellCrop{r},1),1),1:size(cellCrop{r},1),'b')
                if ~isempty(midline)
                    plot(midline(:,1),midline(:,2),'g');
                    plot(contour(:,1),contour(:,2),'m');
                end
                title(qrData{r})
                drawnow
                hold off
                saveas(gca,[oPath qrData{r} '.jpg']);
                imwrite(vec,[oPath qrData{r} '{imageName_straightColor}.jpg']);
                imwrite(vecM,[oPath qrData{r} '{imageName_straightMask}.jpg']);
            else
                fprintf(['Empty Cell.\n']);
                imshow(cellCrop{r},[]);
                title('Empty');
            end
        end
    catch ME
        getReport(ME)
    end
end
%}
%%

%% construct GMM
NC = 3;
options = statset('Display','iter','MaxIter',400);
gmmCarrot = fitgmdist(rootStack,NC,'Options',options);


