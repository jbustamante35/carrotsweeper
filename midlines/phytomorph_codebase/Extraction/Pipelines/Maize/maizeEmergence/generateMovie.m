function [MOV] = generateMovie(imageStack,mask,popN,sidx)
   
   
    for e = 1:size(popN,2)
        tmp = find(popN(:,e));
        popN(:,e) = 0;
        if ~isempty(tmp)
            popN(tmp(1),e) = 1;
        end
    end
    mask = imresize(mask,.25);    
    mask = logical(mask);
    R = regionprops(mask,'Centroid');
    R(sidx) = R;
    
    POPT = logical(zeros(size(mask)));
    
    for e = 2:1:(numel(imageStack)-1)
        I = imread(imageStack{e});
        I = imresize(I,.25);
        fidx = find(popN(e,:));
        POPTN = logical(zeros(size(mask)));
        for f = 1:numel(fidx)
            tmp = imfill(~mask,fliplr((round(R(fidx(f)).Centroid))));
            POPT = POPT | (tmp & mask);
            tmp = imfill(~mask,fliplr(round(R(fidx(f)).Centroid)));
            POPTN = POPTN | (tmp & mask);
            f
        end
        
        out = flattenMaskOverlay(I,mask,.3,'r');
        
        out = flattenMaskOverlay(out,POPT,.4,'b');
        out = flattenMaskOverlay(out,POPTN,.8,'g');
        
        imshow(out,[]);
        title(num2str(e))
        drawnow
        MOV(:,:,:,e-1) = out;
    end
end