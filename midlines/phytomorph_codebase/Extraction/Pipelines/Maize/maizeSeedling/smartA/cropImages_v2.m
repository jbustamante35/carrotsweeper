function [returnI,boundingBoxes] = cropImages_v2(FileList,tester,tester2,PARA)
    fprintf(['..................................\n']);
    fprintf(['..................................\n']);
    
    
    
    I = double(imread(FileList{1}));
    [I,angle] = rectifyImage(I/255);
    I = I*255;
    
    scaleR = (size(I,1)/3280)^-1;
    
    importFunc = @(X)imresize(255*rectifyImage(double(imread(X))/255),[3280 size(I,2)]);
    % create temp data source
    tmpSource = dataSource(FileList,importFunc,0);
    % set the functions to have the temp data sources
    tester.dataSource.dataSource.dataSource = tmpSource;
    tester2.dataSource.dataSource.sourceNozzles{1}.dataSource = tmpSource;
    % create nozzle pointer
    tester.nozzlePtr = nozzlePtr(numel(FileList),1);
    tester3 = dataNozzle(@(X,e0,e1)func_thumbNail(X,[1 4992],false,false),tester,1);
    tester2.dataSource.dataSource.sourceNozzles{2}.dataSource = tester3;
    tester2.nozzlePtr = nozzlePtr(numel(FileList)*3,1);
    tester2.dataSource.dataSource.readsPerSource = 3;
    tester2.dataSource.dataSource.nozzlePtr = nozzlePtr(1,3);
    
    %{
    tester2.dataSource.dataSource.sourceNozzles{1}.dataSource.inMEM = 0;
    tester2.dataSource.dataSource.sourceNozzles{2}.dataSource.inMEM = 0;
    %}

    tester.resetPtr();
    tester2.resetPtr();

    vertStrip = tester.next();
    tester.resetPtr();
    vertStrip = imresize(vertStrip,[1 size(I,2)]);
    R = regionprops(logical(vertStrip),'PixelIdxList','BoundingBox');



    cnt = 1;
    MSK = zeros(size(I,1),size(I,2));
    %while tester2.hasNext()
    for e = 1:numel(R)
        z1 = zeros(size(vertStrip));
        z1(R(cnt).PixelIdxList) = 1;
        %tester.resetPtr();
        horiStrip = tester2.next();
        horiStrip = imresize(horiStrip',[size(I,1) 1],'nearest');
        
        MSK = MSK + horiStrip*z1;
        cnt = cnt + 1;
    end
    
    MSK(1:PARA{1},:) = 0;
    
    
    RT = regionprops(logical(MSK),'PixelIdxList','BoundingBox');
    for e = 1:numel(R)
        %RT(e).BoundingBox(2) = RT(e).BoundingBox(2)*scaleR;
        %RT(e).BoundingBox(4) = RT(e).BoundingBox(4)*scaleR;
        
        RT(e).BoundingBox(1) = RT(e).BoundingBox(1) - PARA{2};
        RT(e).BoundingBox(1) = max([1 RT(e).BoundingBox(1)]);
        RT(e).BoundingBox(3) = RT(e).BoundingBox(3) + 2*PARA{2};
        if (RT(e).BoundingBox(1) + RT(e).BoundingBox(3)) > size(I,2)
            RT(e).BoundingBox(3) = -(RT(e).BoundingBox(1) - size(I,2));
        end
        returnI{e} = imcrop(I,RT(e).BoundingBox);
        boundingBoxes{e} = RT(e).BoundingBox;
    end
    
     
       
    
    fprintf(['..................................\n']);
    fprintf(['..................................\n']);
    
end