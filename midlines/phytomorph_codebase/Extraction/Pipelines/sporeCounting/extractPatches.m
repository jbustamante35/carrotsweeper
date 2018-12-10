function [patchStructure,ntmp] = extractPatches(fileName,patchSZ,areaTHRESH)
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    fprintf(['Start patch extration\n']);tm = clock;
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    
    patchStructure = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find species name from file name
    speciesIDX = strfind(fileName,'sol');
    if isempty(speciesIDX)
        speciesIDX = strfind(fileName,'alt');
    end
    species = fileName(speciesIDX(1):speciesIDX(1)+2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read the image
    fprintf(['start: Reading image file.\n']);tic;
    % read the file
    tmp = imread(fileName);
    fprintf(['end: Reading image file:' num2str(toc) '\n']);
    fprintf(['start: Converting file.\n']);tic;
    % convert to gray scale
    tmp = rgb2gray(tmp);
    % convert to double
    tmp = double(tmp)/255;
    fprintf(['end:  Converting file:' num2str(toc) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % remove background
    fprintf(['start: Removing background.\n']);tic;
    % resize for background calculation
    itmp = imresize(tmp,.25);
    % filter the background
    BK = imfilter(itmp,fspecial('disk',101),'replicate');
    % upsize the backgroud
    BK = imresize(BK,size(tmp));
    % remove the background
    ntmp = tmp - BK;
    % normalize the image
    ntmp = bindVec(ntmp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % perform adaptiave threshold calculation
    T = adaptthresh(ntmp,.4,'ForegroundPolarity','dark');
    % perform calculation
    msk = ~imbinarize(ntmp,T);
    % remove small objects
    msk = bwareaopen(msk,areaTHRESH);
    fprintf(['end:  Removing background:' num2str(toc) '\n']);
    fprintf(['start: Finding and clipping objects.\n']);tic;
    % find regions
    R = regionprops(msk,'Area','PixelIdxList','Centroid','Eccentricity');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % multi-level count
    % store image patches in "database" for later use
    cnt = 1;
    for iter = 1:4
        % count the objects based on area
        fidx = count([R.Area]);
        fidx1 = count([R.Eccentricity]);
        % find the common objects
        sR = R(fidx==1 & fidx1 == 1);
        % remove from the set
        R(fidx==1 & fidx1 == 1) = [];
        
        
        for o = 1:numel(sR)
            % round the centroid
            R(o).Centroid = round(sR(o).Centroid);
            try
                patchStructure(cnt).patch = ntmp(sR(o).Centroid(2)-patchSZ:sR(o).Centroid(2)+patchSZ,sR(o).Centroid(1)-patchSZ:sR(o).Centroid(1)+patchSZ);
                patchStructure(cnt).mask = msk(sR(o).Centroid(2)-patchSZ:sR(o).Centroid(2)+patchSZ,sR(o).Centroid(1)-patchSZ:sR(o).Centroid(1)+patchSZ);
                patchStructure(cnt).patchLocation = sR(o).Centroid;
                patchStructure(cnt).species = species;
                cnt = cnt + 1;
            catch
                
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    fprintf(['end:Finding and clipping objects:' num2str(toc) '\n']);
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    fprintf(['End patch extration.' num2str(etime(clock,tm)) '\n']);
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
end