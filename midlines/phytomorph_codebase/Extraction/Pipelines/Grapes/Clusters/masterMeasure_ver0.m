function [] = masterMeasure_ver0(FileName,clusterForest,oPath)

    [p,nm,ext] = fileparts(FileName);

    oI = double(imread(FileName));
    
    % major mask
    clusterTree_1 = clusterForest.getTree(3);
    k1 = clusterTree_1.clusterImage(FileName);
    Mask = k1 == 3;
    [phDocument] = grapeClusterMetric_ver1(Mask);
    displayPhenotypes(oI,Mask,phDocument,oPath,nm);
    
    
    
    
    
    
    % subcluster
    clusterTree_2 = clusterForest.getTree(1);
    k2 = clusterTree_2.clusterImage(FileName);
    clusterVec = [9 12 13 15 16];
    for v = 1:numel(clusterVec)
        clusterMask(:,:,v) = k2 == clusterVec(v);
        [phDocument2{v}] = grapeClusterMetric_ver1(clusterMask(:,:,v));
        displayPhenotypes(oI,Mask,phDocument2{v},oPath,nm);
    end
    CL = {'r' 'g' 'b' 'y' 'm'};
    out = oI/255;
    for v = 1:numel(clusterVec)
        out = flattenMaskOverlay(out,clusterMask(:,:,v),.2,CL{v});
    end
    image(out);
    axis off
    hold on
    for e = 1:numel(phDocument2)
        plot(phDocument2{e}.massDistributionVer,1:numel(phDocument2{e}.massDistributionVer),CL{e});
        plot(1:numel(phDocument2{e}.massDistributionHor),phDocument2{e}.massDistributionHor,CL{e});
    end
    
    
    
    
    
    
  %{
    Mask = 
    
    k = forest.clusterImage(FileList{i});
    for e = 1:size(k,3)
        rgb = label2rgb(k(:,:,e));
        imshow(cat(2,oI/255,double(rgb)/255),[]);
        drawnow
        waitforbuttonpress
    end
    %}
end