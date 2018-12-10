function [M,out] = label2Mask(I,L,maskV,CL,disp)
    out = I;
    for e = 1:numel(maskV)
        M(:,:,e) = zeros(size(I,1),size(I,2));
        for l = 1:numel(maskV{e})
            M(:,:,e) = logical(M(:,:,e)) | L == maskV{e}(l);
        end
        if disp
            out = flattenMaskOverlay(out,logical(M(:,:,e)),.4,CL{e});
        end
    end
    if disp
        imshow(out,[]);
    end
end