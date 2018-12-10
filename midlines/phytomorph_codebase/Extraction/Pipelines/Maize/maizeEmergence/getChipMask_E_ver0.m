function [M ampSig] = getChipMask_E_ver0(stack)
    disp = 0;
    M = ((stack(:,:,1,:) < .08) | (stack(:,:,1,:) > .8)) & stack(:,:,2,:) > .25;
    for e = 1:size(M,4)
        M(:,:,1,e) = imfill(M(:,:,1,e),'holes');
        %M(:,:,1,e) = imclearborder(M(:,:,1,e));
        vec = stack(:,:,2,e).*M(:,:,1,e);
        ampSig(e) = std(vec(:));
        if disp
            imshow(M(:,:,1,e),[])
            drawnow
        end
    end
    
end