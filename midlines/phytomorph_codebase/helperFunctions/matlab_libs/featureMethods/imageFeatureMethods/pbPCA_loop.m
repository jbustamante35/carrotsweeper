function [Y Tout] = pbPCA_loop(patchStack,domain)
    disp = 0;
    num_tr = size(patchStack,ndims(patchStack));
    patchSize = size(patchStack);
    patchSize = patchSize(1:end-1);
    for e = 1:num_tr
        T = [eye(2) ((patchSize-1)/2)'];
        [Tout(:,:,e) Y(:,e)] = pbPCA(patchStack(:,:,e),T,domain,@gt);
        if disp
            imshow(patchStack(:,:,e),[]);
            hold on
            quiver(Tout(2,end),Tout(1,end),Tout(2,1),Tout(1,1));
        end
    end
end