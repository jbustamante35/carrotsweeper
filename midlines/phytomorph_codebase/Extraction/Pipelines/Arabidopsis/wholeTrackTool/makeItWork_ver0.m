function [VEC] = makeItWork_ver0(lFlow,I,re,sz)
    lFlow = bsxfun(@minus,lFlow,lFlow(:,:,:,1));

    SNIP = 100;
    relFlow = [];
    for e = 1:size(lFlow,4)
        for d = 1:size(lFlow,3)
            relFlow(:,:,d,e) = imresize(lFlow(SNIP:(end-(SNIP-1)),SNIP:(end-(SNIP-1)),d,e),re);
        end
    end
    
    F1 = [];
    for e = 1:size(relFlow,4)
        tmp = [];
        for d = 1:size(relFlow,3)
            tmp = [tmp;im2col(relFlow(:,:,d,e),[sz sz],'sliding')];
        end
        F1 = [F1;tmp];
        e
    end
    
    reI = imresize(I(SNIP:(end-(SNIP-1)),SNIP:(end-(SNIP-1))),re);
    F2 = im2col(reI,[sz sz],'sliding');
    VEC = [F1;F2];
end