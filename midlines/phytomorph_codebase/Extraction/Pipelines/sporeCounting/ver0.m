FilePath = '/mnt/tetra/nate/Spore/';
FilePath = '/mnt/tetra/nate/Spore/subset4counting/';
FileList = {};
FileExt = {'JPG'};
FileList = gdig(FilePath,FileList,FileExt,1);
%%
FileList = FileList(randperm(numel(FileList)));
%%
searchPath = '/iplant/home/nmiller/subset4counting%';
searchPath = '/iplant/home/nmiller/sporeMaster/All4Nathan%';
CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' searchPath '''"']; 
[o,r] = system(CMD);
[r] = parseRecords(r);
FileList = {};
for e = 1:numel(r)
    FileList{e} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
end
%%

%%
patchSZ = 60;
oPath = '/mnt/tetra/nate/Spore/returns/return2/';
mkdir(oPath);
tic
masterTable = table;
for e = 1:numel(FileList)
    fprintf(['**************************************************************************************\n']);
    fprintf(['**************************************************************************************\n']);
    fprintf(['Start image analysis\n']);tm = clock;
    fprintf(['**************************************************************************************\n']);
    fprintf(['**************************************************************************************\n']);
    [pth,fileName,ext] = fileparts(FileList{e});
    
    %[out,pTable] = masterSpore(FileList{e},patchSZ,'');
    
    [patchStructure,nI] = extractPatches(FileList{e},patchSZ,100);
    
    [patchStructure,pTable] = extractPixelListsFromPatches(patchStructure,size(nI),FileList{e});
    
    masterTable = [masterTable;pTable];
    writetable(pTable,[oPath fileName '___phenotypeTable.txt'])
    [out] = paintImage(nI,patchStructure,oPath,fileName);
    
    fprintf(['**************************************************************************************\n']);
    fprintf(['**************************************************************************************\n']);
    fprintf(['End image analysis.' num2str(etime(clock,tm)) '\n']);
    fprintf(['**************************************************************************************\n']);
    fprintf(['**************************************************************************************\n']);
end
toc
%%
[FileList] = issueBulkTicket(FileList);
%%
FileList(1) = [];
%%
FileList{BAD(2)}
%%
func = cFlow('masterSpore');
func.setMCRversion('v930');
patchSZ = [60];
imageOut = {};
tableOut = {};
for e = 1:numel(FileList)
    [imageOut{e},tableOut{e}] = func(FileList{e},patchSZ,'');
    e
end
auth = readtext('/mnt/spaldingdata/nate/auth.iplant');
auth = auth{1};
func.submitDag(auth,150,150);



%%  load master table
masterTable = table;
iPath = '/mnt/tetra/nate/sporeReturn/';
mkdir(iPath);
tm = [];
BAD = [];
T = {};
for e = 1:numel(tableOut)
    tic
    try
        T{e} = cFlowLoader(tableOut{e});
        writetable(T{e},[iPath num2str(e) '.csv']);
        %masterTable = [masterTable;T];
        tm(e) = toc;
        %mean(tm)*(numel(tableOut)-e)
        fprintf(['Done loading:' num2str(e) '\n']);
        BAD(e) = 0;
    catch ME
        ME
        BAD(e) = 1;
    end
end
%%
%%  try load after crash
FilePath = '/mnt/tetra/nate/sporeReturn/';
csvFileList = {};
FileExt = {'csv'};
csvFileList = gdig(FilePath,csvFileList,FileExt,1);
%% 
%SZ = zeros(numel(csvFileList),1);
str = 1;
S = zeros(sum(SZ)-numel(csvFileList),39);
%%
close all
cnt = 1;
for e = 1:numel(csvFileList)
    
    tic
    data = readtext(csvFileList{e});
    try
        %
        NM1{cnt} = data(:,1);
        NM2{cnt} = data(:,7);



        data(1,:) = [];
        data(:,1) = [];
        data(:,7) = [];



        data = cell2mat(data);
        stp = str + size(data,1) - 1;
        S(str:stp,:) = data;

        str = stp + 1;
        %SZ(e) = size(data,1);
        toc
        e
        delta(e) = toc;
        tm(e) = ((numel(csvFileList) - e)*mean(delta))/60/60;
        %SZ(e) = size(data,1);
        plot(tm)
        drawnow
        cnt = cnt + 1;
    catch
        
    end
end
%%
tot = 0;
for e = 1:numel(NM1)
    
    tot = tot + size(NM2{e},1);
    e
    
end
%%
mtf = zeros(numel(csvFileList),1);
for e = 1:numel(csvFileList)
    
    tic
    data = readtext(csvFileList{e});
    
    H = data(:,1);
    

    data(1,:) = [];
    data(:,1) = [];
    data(:,7) = [];



    data = cell2mat(data);
    
    if size(data,1) == size(H,1)-1
        mtf(e) = 1;
    else
        mtf(e) = 0;
        break
    end
    
    toc
end
%% split into genotypes
fidx1 = strcmp(masterTable.species,'alt');
fidx2 = strcmp(masterTable.species,'sol');
altTable = masterTable(fidx1,:);
solTable = masterTable(fidx2,:);
%% look at distributions
sam = 10;
plot3(features(1:sam:end,1),features(1:sam:end,2),features(1:sam:end,3),'.');
%% stack vars for alt table
toUse = altTable.Properties.VariableNames;
toUse([1 6 7 8]) = [];
IDX = setdiff(1:numel(altTable.Properties.VariableNames),[1 6 7 8]);

dataVec = double(table2array(altTable(:,IDX)));
fName1 = altTable.fileName;
aIDX = [altTable.columnIndex altTable.rowIndex];
rmidx = find(any(isnan(dataVec),2));
fName1(rmidx) = [];
dataVec(rmidx,:) = [];
aIDX(rmidx,:) = [];

UQ = unique(fName1);
[S C U E L ERR LAM] = PCA_FIT_FULL(dataVec,size(dataVec,2));
%% stack vars for sol table
toUse = solTable.Properties.VariableNames;
toUse([1 6 7 8]) = [];
IDX = setdiff(1:numel(solTable.Properties.VariableNames),[1 6 7 8]);

dataVec = double(table2array(solTable(:,IDX)));
fName1 = solTable.fileName;
aIDX = [solTable.columnIndex solTable.rowIndex];
rmidx = find(any(isnan(dataVec),2));
fName1(rmidx) = [];
dataVec(rmidx,:) = [];
aIDX(rmidx,:) = [];

UQ = unique(fName1);
[S C U E L ERR LAM] = PCA_FIT_FULL(dataVec,size(dataVec,2));
%% 
oPath = '/mnt/tetra/nate/retSporeTable/';
mkdir(oPath)
writetable(altTable,[oPath 'altTable.csv']);
writetable(countTable,[oPath 'count_altTable.csv']);
%% get clusters
oPath = '/mnt/tetra/nate/retSporeSOL/';
mkdir(oPath);
TOTU = numel(UQ);
%TOTU = 20;
REP = 10;
countTable = table;
for clusterK = 4%3:5%4%3:5
    for fNum = 4%2:5%4%2:5

        options = statset('Display','iter','MaxIter',300);
        GMModel = fitgmdist(C(:,1:fNum),clusterK,'Options',options,'RegularizationValue',0.00001,'Replicates',REP);

        kidx = GMModel.cluster(C(:,1:fNum));
        %
        CL = {'r.' 'g.' 'b.' 'c.' 'm.'};
        
        for u = 1:TOTU
            fidx = find(strcmp(fName1,UQ{u}));
            subK = kidx(fidx);
            
            for k = 1:clusterK
                N = sum(subK == k);
                if k == 1
                    ptr = size(countTable,1)+1;
                end
                countTable{ptr,'fileName'} = UQ(u);
                countTable{ptr,['group' num2str(k)]} = N;
            end
            countTable{ptr,'TOTAL'} = numel(subK);
        end
        
        %{
        parfor u = 1:TOTU
            close all
            fidx = find(strcmp(fName1,UQ{u}));
            I = imread(UQ{u});
            subK = kidx(fidx);
            subIDX = aIDX(fidx,:);
            out = I;
            
            %imshow(I,[]);
            %hold on
            for k = 1:numel(CL)
                skidx = find(subK == k);
                msk = zeros(size(I,1),size(I,2));
                tIDX = sub2ind(size(msk),round(subIDX(skidx,2)),round(subIDX(skidx,1)));
                msk(tIDX) = 1;
                msk = imdilate(msk,strel('disk',2,0));
                out = flattenMaskOverlay(out,logical(msk),.9,CL{k}(1));
                
                plot(subIDX(skidx,1),subIDX(skidx,2),CL{k});
            end
            [~,nm] = fileparts(UQ{u});
            %imshow(out,[]);
            oPathT = [oPath 'C_' num2str(clusterK)];
            mkdir(oPathT);
            oPathT = [oPath 'C_' num2str(clusterK) filesep 'F_' num2str(fNum) filesep];
            mkdir(oPathT);
            oName = [oPathT nm '.jpg'];
            imwrite(out,oName);
            oName
        end
        %}
        
        
    end
end
%% cluster tables
features = [solTable.meanHyphaDistance,solTable.stdHyphaDistance,solTable.sporeEccentricity];
GMModel = fitgmdist(features,3);
%% cluster on features
features = [masterTable.meanHyphaDistance,masterTable.stdHyphaDistance,masterTable.sporeEccentricity];
GMModel = fitgmdist(features,3);
%% display clusters
kIDX = cluster(GMModel,features);
close all
gscatter(features(:,1),features(:,2),kIDX)
hold on
fcontour(@(x1,x2)pdf(GMModel,[x1 x2]),[min(features(:,1)) max(features(:,1)) min(features(:,2)) max(features(:,2))])
%% label the "alt" species as 0 on the S vector
UQ = unique(masterTable.species);
for u = 1:numel(UQ)
    speciesGroup = masterTable.species;
    speciesGroup = strcmp(speciesGroup,UQ{u});
    IDX = find(speciesGroup==1);
    ksdensity(masterTable.totalHypaArea(IDX));
    waitforbuttonpress
end
%{
%%
    %kidx = kmeans(masterTable.totalSporeArea(solIDX),2);
    % find the sol greater than 6000 - list those as 1 or 2
    kidx = (masterTable.totalSporeArea(solIDX) > 600) + 1;
    % create master kidx
    kIDX = zeros(size(masterTable,1),1);
    % create a feature vector for the SOL speciea
    solfeatureVector = [masterTable.meanHyphaDistance(solIDX(kidx==2)),masterTable.stdHyphaDistance(solIDX(kidx==2))];
    % fit two clusters to the SOL
    GMModel = fitgmdist(solfeatureVector,2);
    % cluster the SOL into the two groups
    SOLDA = cluster(GMModel,solfeatureVector);
    %SOLDA = kmeans(solfeatureVector,2);
    % label the SOL as 1 or 2
    kIDX(solIDX(kidx==1)) = 1;
    kIDX(solIDX(kidx==2)) = SOLDA+1;
end
%}
%% find the sol
solIDX = find(S==1);
%kidx = kmeans(masterTable.totalSporeArea(solIDX),2);
% find the sol greater than 6000 - list those as 1 or 2
kidx = (masterTable.totalSporeArea(solIDX) > 600) + 1;
% create master kidx
kIDX = zeros(size(masterTable,1),1);
% create a feature vector for the SOL speciea
solfeatureVector = [masterTable.meanHyphaDistance(solIDX(kidx==2)),masterTable.stdHyphaDistance(solIDX(kidx==2))];
% fit two clusters to the SOL
GMModel = fitgmdist(solfeatureVector,2);
% cluster the SOL into the two groups
SOLDA = cluster(GMModel,solfeatureVector);
%SOLDA = kmeans(solfeatureVector,2);
% label the SOL as 1 or 2
kIDX(solIDX(kidx==1)) = 1;
kIDX(solIDX(kidx==2)) = SOLDA+1;
%% 
UQ = unique(kIDX);
UQ(isnan(UQ)) = [];
close all
zoomLevel = 250;
for u = 3:4%1:numel(UQ)
    fidx = find(kIDX==UQ(u));
    %fidx = fidx(randperms(numel(fidx)));
    tmpName = '';
    oldName = '';
    for f = 1:100
        try
            tmpName = masterTable(fidx(f),'fileName');
            tmpName = tmpName.fileName{1};


            sporearea = masterTable(fidx(f),'totalSporeArea');
            sporearea = sporearea.totalSporeArea;
            ypos = masterTable(fidx(f),'rowIndex');
            xpos = masterTable(fidx(f),'columnIndex');
            xpos = xpos.columnIndex;
            ypos = ypos.rowIndex;

            if ~strcmp(tmpName,oldName)
                tmpI = imread(tmpName);
                imshow(tmpI,[]);
            end

            hold on
            plot(xpos,ypos,'c*');
            axis([xpos - zoomLevel xpos + zoomLevel ypos - zoomLevel ypos + zoomLevel]);
            gidx0 = strfind(tmpName,'EXP');
            gidx1 = strfind(tmpName,'/');

            speciesIDX = strfind(tmpName,'sol');
            if isempty(speciesIDX)
                speciesIDX = strfind(tmpName,'alt');
            end

            species = tmpName(speciesIDX(1):speciesIDX(1)+2);

            gidx1(gidx1 < gidx0(1)) = [];
            folderName = tmpName(gidx0(1):(gidx1(1)-1));
            title(['species:' species '--folder:' folderName  '---group:' num2str(u) '--spore:' num2str(f) '--area:' num2str(sporearea)]);
            waitforbuttonpress
            oldName = tmpName;
        end
       
    end
    
end
%%
close all
pSZ = [300 300];
for e = 1:numel(FileList)
    tmp = imread(FileList{e});
    tmp = rgb2gray(tmp);
    tmp = double(tmp)/255;
    itmp = imresize(tmp,.25);
    BK = imfilter(itmp,fspecial('disk',101),'replicate');
    BK = imresize(BK,size(tmp));
    ntmp = tmp - BK;
    ntmp = bindVec(ntmp);
    msk = ntmp < graythresh(ntmp);
    T = adaptthresh(ntmp,.4,'ForegroundPolarity','dark');
    msk = ~imbinarize(ntmp,T);
    msk = bwareaopen(msk,50);
    R = regionprops(msk,'Area','PixelIdxList','Centroid');
    fidx = count([R.Area]);
    R = R(fidx==1);
    nmsk = zeros(size(msk));
    for o = 1:numel(R)
        nmsk(R(o).PixelIdxList) = 1;
    end
    
    out = flattenMaskOverlay(ntmp,msk);
    out = flattenMaskOverlay(out,logical(nmsk),.8,'g');
    %G(:,:,e) = tmp;
    s = [[bindVec(tmp(1:pSZ(1),1:pSZ(2))) bindVec(ntmp(1:pSZ(1),1:pSZ(2)))];...
    [bindVec(ntmp(1:pSZ(1),1:pSZ(2))) bindVec(out(1:pSZ(1),1:pSZ(2)))]];
    imshow(out,[]);
    drawnow
end
%%

close all
pSZ = [300 300];
patchSZ = [60];
cnt = 1;
disp = false;
pS = [];
pM = [];
imageNumber = [];
patchLocation=  [];
G = [];
for e = 1:100%numel(FileList)
    tmp = imread(FileList{e});
    tmp = rgb2gray(tmp);
    tmp = double(tmp)/255;
    itmp = imresize(tmp,.25);
    BK = imfilter(itmp,fspecial('disk',101),'replicate');
    BK = imresize(BK,size(tmp));
    ntmp = tmp - BK;
    ntmp = bindVec(ntmp);
    T = adaptthresh(ntmp,.4,'ForegroundPolarity','dark');
    msk = ~imbinarize(ntmp,T);
    msk = bwareaopen(msk,50);
    R = regionprops(msk,'Area','PixelIdxList','Centroid');
    fidx = count([R.Area]);
    R = R(fidx==1);
    nmsk = zeros(size(msk));
    G(:,:,e) = ntmp;
    for o = 1:numel(R)
        nmsk(R(o).PixelIdxList) = 1;
        R(o).Centroid = round(R(o).Centroid);
        try
            pS(:,:,cnt) = ntmp(R(o).Centroid(2)-patchSZ:R(o).Centroid(2)+patchSZ,R(o).Centroid(1)-patchSZ:R(o).Centroid(1)+patchSZ);
            pM(:,:,cnt) = msk(R(o).Centroid(2)-patchSZ:R(o).Centroid(2)+patchSZ,R(o).Centroid(1)-patchSZ:R(o).Centroid(1)+patchSZ);
            imageNumber(cnt) = e;
            patchLocation(cnt,:) = R(o).Centroid;
            cnt = cnt + 1;
        catch
        end
        cnt
        e
    end
    if disp
        out = flattenMaskOverlay(ntmp,msk);
        out = flattenMaskOverlay(out,logical(nmsk),.8,'g');
        imshow(out,[]);
        drawnow
    end
end
%%

%%
close all
cp = (size(pM)+1)/2;
masterMask = zeros(size(G));
masterMask2 = zeros(size(G));
for e = 1:size(pS,3)
    try
        
        tmpMM = masterMask(:,:,imageNumber(e));
        tmpMMM = masterMask2(:,:,imageNumber(e));
        
        
        tmpLoc = flipdim(patchLocation(e,:),2) - patchSZ;
        para.scales.value = 2;
        para.resize.value = 1;
        tmpM = imfill(pM(:,:,e),'holes');
        msk = imfill(~tmpM,cp(1:2)) == 1 & pM(:,:,e) == 1;
        [K] = surKur(pS(:,:,e),para);
        K = bindVec(K(:,:,1));

        midx = find(msk);
        
        what = [];
        [what(:,1) what(:,2)] = ind2sub(size(tmpM),midx);
        what = bsxfun(@plus,what,tmpLoc);
        what = sub2ind([size(G,1),size(G,2)],what(:,1),what(:,2));
        tmpMMM(what) = 1;
        
        
        K(midx) = mean(K(:));
        d = K > graythresh(K);

        T = adaptthresh(K,.4);
        d = imbinarize(K,T);


        Rh = regionprops(d,'PixelIdxList');
        Rs = regionprops(imdilate(msk,strel('disk',3,0)),'PixelIdxList');
        H = zeros(size(d));
        for h = 1:numel(Rh)
            if ~isempty(intersect(Rs(1).PixelIdxList,Rh(h).PixelIdxList))
                H(Rh(h).PixelIdxList) = 1;
                
                what = [];
                [what(:,1) what(:,2)] = ind2sub(size(tmpM),Rh(h).PixelIdxList);
                what = bsxfun(@plus,what,tmpLoc);
                what = sub2ind([size(G,1),size(G,2)],what(:,1),what(:,2));
                tmpMM(what) = 1;
            end
        end
        H = logical(H);

        msk1 = d == 1 & msk == 1;
        msk2 = d == 1 & msk == 0;
        msk3 = d == 0 & msk == 1;
        out = flattenMaskOverlay(pS(:,:,e),msk1,.2,'r');
        out = flattenMaskOverlay(out,msk2,.2,'g');
        out = flattenMaskOverlay(out,msk3,.2,'b');
        out = flattenMaskOverlay(out,H,.8,'m');
        imshow(out,[])
        drawnow
        %pause(.5);
    catch
    end
    masterMask(:,:,imageNumber(e)) = tmpMM;
    masterMask2(:,:,imageNumber(e)) = tmpMMM;
  
end
%%
oPath = '/mnt/tetra/nate/Alt.sea/return/colored/';
for e = 1:size(masterMask,3)
    out = flattenMaskOverlay(G(:,:,e),logical(masterMask(:,:,e)),.3,'g');
    out = flattenMaskOverlay(out,logical(masterMask2(:,:,e)),.3,'b');
    imwrite(out,[oPath num2str(e) '.jpg']);
    imshow(out,[]);
    drawnow
end
%%
oPath = '/mnt/tetra/nate/Alt.sea/return/raw/';
for e = 1:size(masterMask,3)
    imwrite(G(:,:,e),[oPath num2str(e) '.jpg']);
end