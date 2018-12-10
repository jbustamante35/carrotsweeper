% scan for image files
FilePath = '//mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/photos/';
FilePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/photos/growthchamber12.7.15/';
FilePath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/seedlingData/growthchamber12.7.15/';

FileList = {};
FileExt = {'jpeg','tif'};
verbose = 1;
FileList = gdig(FilePath,FileList,FileExt,verbose);
%% view images
for e = 1:numel(FileList)
    I = imread(FileList{e});
    imshow(I,[]);
    drawnow
end
%% loop over images
oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/return/seedlings/output/';
for e = 1:numel(FileList)
    singleSeedlingImage(FileList{e},50,5,100,100,4,oPath);
end
%% find the top container line and pour data into STORE
STORE = {};
EXT = 100;
threshSIG = 5;
for e = 1:numel(FileList)
    % read the image
    I = imread(FileList{e});
    % make gray scale
    G = rgb2gray(I);
    G = imfilter(G,fspecial('gaussian',[13 13],4));
    % find edge
    E = edge(G);
    % integrate
    sig = sum(E,1);
    % smooth
    sig = imfilter(sig,fspecial('average',[1 50]),'replicate');
    BLOCK = sig < threshSIG;
    BLOCK = bwareaopen(BLOCK,50);
    BLOCK = imclose(BLOCK,strel('disk',100));
    eBLOCK = imerode(BLOCK,strel('disk',[EXT]));
    MASK = repmat(eBLOCK,[size(I,1) 1]);
    R = regionprops(~MASK,'BoundingBox');
    
    
    imshow(I,[]);
    hold on
    plot(double(BLOCK)*1000,'r')
    plot(double(eBLOCK)*1000,'g');
    hold off
    drawnow
    
    close all
    for e = 1:numel(R)
        tmpD = imcrop(I,R(e).BoundingBox);
        tmpD(round(1:.45*size(I,1)),:,:) = [];
        imshow(tmpD,[])
        drawnow
        pause(1);
        STORE{end+1} = tmpD;
    end
end
%% analysis of STORE and subCROP the STORE
close all
eT = 120;
OFFSET = 80;
sigFILL = 1100;
cSTORE = {};
for e = 1:numel(STORE)
    G = rgb2gray(STORE{e});
    G = imfilter(G,fspecial('disk',15),'replicate');
    SZ = size(STORE{e});
    E = edge(G);
    sig = sum(E,2);
    sig = sig > eT;
    sig(1:sigFILL) = 0;
    %[J,nidx] = max(sig);
    nidx = find(sig);
    if ~isempty(nidx)
        nidx = nidx(1);
        cSTORE{end+1} = STORE{e}(1:(nidx-OFFSET),:,:);
        imshow(STORE{e});
        hold on
        plot(1:SZ(2),nidx*ones(1,SZ(2)),'r')
        plot(1:SZ(2),(nidx-OFFSET)*ones(1,SZ(2)),'g')
        hold off
        drawnow
    end
end
%% view cSTORE
oPath = '/mnt/spaldingdata/nate/mirror_images/maizeData/hirsc213/photos/returnTemp/';
close all
SNIP = 5;
thresP = .3;
for e = 1:numel(cSTORE)
    tmpD = cSTORE{e};
    SZ = size(tmpD);
    MASK = getMASK_ver0(tmpD);
    if sum(MASK(:))/prod(size(MASK)) < thresP
        tmpMASK = padarray(MASK, [300 0], 'replicate', 'post');
        SKEL = bwmorph(tmpMASK,'thin',inf);
        SKEL = SKEL(1:size(SKEL,1),:);

        [r c] = find(SKEL);
        %{
        R = regionprops(MASK,'Area','PixelIdxList');
        [J midx] = max([R.Area]);
        MASK = zeros(size(MASK));
        MASK(R(midx).PixelIdxList) = 1;
        MASK = logical(MASK);
        %}

        sig = sum(MASK,2);
        fidx = find(sig);
        HEIGHT = fidx(1);
        dBIOMASS = sum(MASK(:));

        EP = imfilter(double(SKEL),ones(3,3));
        [re ce] = find(EP == 2 & SKEL);

        baseMASK = sum(MASK(end-SNIP:end,:),1);
        basePoint(1) = mean(find(baseMASK));
        basePoint(2) = size(MASK,1);
        % find skeleton
        [x y] = find(SKEL);
        DP = [x y]';
        T = Radjacency(DP,3);


        % single trace
        pathcost = [];
        path = {};
        [idx(1)] = snapTo(DP',[fliplr(basePoint)]);
        for i = 1:numel(re)
            [idx(2)] = snapTo(DP',[re(i) ce(i)]);
            [path{i} , pathcost(i)]  = dijkstra(T , idx(1) , idx(2));
            %{
            imshow(MASK,[]);
            hold on
            plot(DP(2,path{i}),DP(1,path{i}),'r')
            plot(DP(2,idx(1)),DP(1,idx(1)),'r*')
            plot(DP(2,idx(2)),DP(1,idx(2)),'r*')
            hold off
            drawnow
            waitforbuttonpress
            %}
        end
        pathcost(isinf(pathcost)) = 0;
        ridx = find(pathcost==0);
        pathcost(ridx) = [];
        path(ridx) = [];
        [J,midx] = max(pathcost);








        out = flattenMaskOverlay(double(cSTORE{e})/255, MASK,.65,'g');
        imshow(out,[])
        hold on
        plot(1:SZ(2),HEIGHT*ones(1,SZ(2)),'g')
        %plot(c,r,'r.')
        plot(ce,re,'b*')
        for i = 1:numel(pathcost)
            plot(DP(2,path{i}),DP(1,path{i}),'r','LineWidth',3)  
        end
        plot(DP(2,path{midx}),DP(1,path{midx}),'k','LineWidth',3)
        hold off
        drawnow 
    else
        imshow(tmpD,[]);
    end
    saveas(gca,[oPath num2str(e) '.tif']);
end
%%