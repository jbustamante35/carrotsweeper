function [S] = preLoadKernelImages(Stack,CN,BOX,SKIP,RE)
    idxVec = 1:SKIP:numel(Stack);
    parfor e = 1:numel(idxVec)
        tmp = getKernelImage(Stack{e},CN,BOX);
        S(:,:,:,e) = imresize(tmp,RE);
        fprintf(['Done with :' num2str(e) ':' num2str(numel(Stack)/SKIP) '\n']);
    end
end