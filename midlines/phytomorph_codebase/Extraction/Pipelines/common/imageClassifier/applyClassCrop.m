function [] = applyClassCrop(fileName,T)
    I= imread(fileName);
    I = imread(fileName);
    I = imresize(I,[50 50]);
    
    type = find(T.class.predict(I));
    
    % percent to reduce the image size
    reducPer = .1;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read training image stack - X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(['start reading train images \n']);
    [imageVec] = readImagesFromList({fileName},reducPer);
    fprintf(['stop reading train images \n']);

    
    
    [XdataCOL] = stackColumnData(imageVec,[],true);
         

    
    
    for cropBox = 1:numel(T.compress{type})
        
        [XdataCOL_COMP] = compressForCropBoxMethod(XdataCOL,T.compress{type}{cropBox}.columnCompressMean,T.compress{type}{cropBox}.columnCompressFrame);
        
        
        columnLabels = T.drawers{type}{cropBox}.columnNet.predict(XdataCOL_COMP);
       
        [~,columnLabels] = max(columnLabels{1},[],1);
        MSK = mod(columnLabels,2)==0;
        
        MSK = repmat(MSK,[size(imageVec,1) 1]);
        
        
        
        [XdataROW] = stackRowReSizeData(imageVec,MSK,200);
        [XdataROW_COMP] = compressForCropBoxMethod(XdataROW,T.compress{type}{cropBox}.rowCompressMean,T.compress{type}{cropBox}.rowCompressFrame);
        
        rowLabels = T.drawers{type}{cropBox}.rowNet.predict(XdataROW_COMP);
        
        
        for k = 1:numel(rowLabels)
            [~,tmprowLabels] = max(rowLabels{k},[],1);
            tmpMSK = mod(tmprowLabels,2)==0;
            MASKS{k} = double(tmpMSK)'*double((mod(columnLabels,2)==0));
            
            
        end
        
        
    end
end

%{

 	% scan base path
    cdir = dir(rootPath);
    cdir(1:2) = [];
    masterList = {};
    typeList = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the image lists
    for e = 1:numel(cdir)
        FilePath = [rootPath '/' cdir(e).name '/imageStacks/'];
        FileList = {};
        FileExt = {'tif','TIF','nef','NEF','PNG','png','jpeg','JPG','jpg','tiff'};
        FileList = gdig(FilePath,FileList,FileExt,1);
        masterList = cat(2,masterList,FileList);
        typeList = cat(1,typeList,e*ones(numel(FileList),1));
    end

    applyClassCrop(masterList{200},T)

%}