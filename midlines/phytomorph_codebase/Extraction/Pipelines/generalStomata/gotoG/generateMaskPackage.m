function [fullDataOut] = generateMaskPackage(trainFileList,toLoad,dilateValue,basisU,basisE)

   
    toLoad = min(toLoad,numel(trainFileList));


    for e = 1:toLoad
        data = load(trainFileList{e},'T','rI','customData','uix');
        
        data = prepareData(data,basisU,basisE);
        
        
        %{
        tmpM = zeros(size(data.rI));
        tmpM(data.customData{2}) = 1;
        R = regionprops(logical(tmpM),'boundingBox');
        subI = imcrop(data.rI,R(1).BoundingBox);
        %}
        
        subI = data.rI;
        uix = round(data.uix);
        tmpY = zeros(size(subI));
        for pt = 1:size(uix,1)
            tmpY(uix(pt,1),uix(pt,2)) = 1;
        end

        tmpY = imdilate(tmpY,strel('disk',dilateValue,0));
        out = flattenMaskOverlay(subI,logical(tmpY));

        data.mask = logical(tmpY);
        
        imshow(out,[]);
        hold on
        plot(data.uix(:,2),data.uix(:,1),'y*')
        drawnow
        hold off
        
        if e == 1
            fullDataOut.C = zeros([size(data.C) numel(trainFileList)]);
            fullDataOut.M = zeros([size(data.mask) numel(trainFileList)]);
            fullDataOut.T = zeros([size(data.T) numel(trainFileList)]);
            fullDataOut.oI = zeros([size(subI) numel(trainFileList)]);
        end
        
        
        fullDataOut.C(:,:,e) = data.C;
        fullDataOut.M(:,:,e) = data.mask;
        fullDataOut.T(:,:,:,:,e) = data.T;
        fullDataOut.oI(:,:,e) = subI;

        
        e
        
    end
end