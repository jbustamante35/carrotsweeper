function [B] = prep(mv,reSize,gridSize)
    tB = reSizeBlocks(mv(:,:,:,1),reSize,gridSize);
    
    B = zeros([size(tB) size(mv,4)],'single');
    for i = 1:size(mv,4)
        I = mv(:,:,:,i);
        B(:,:,i) = reSizeBlocks(I,reSize,gridSize);
        i
    end
end


function [tBlock] = reSizeBlocks(I,reSize,gridSize)
    tBlock = [];
    for r = 1:numel(reSize)
        rI = imresize(I,reSize(r));
        cBlock = [];
        for k = 1:size(rI,3)
            %fun = @(m)func(m,k);
            %tmp(:,:,k) = colfilt(rI(:,:,k),gridSize,'sliding',fun);
            cBlock = cat(3,cBlock,im2colF(rI(:,:,k),gridSize(1:2),[1 1]));
            %tmp(:,k) = fun(block);
            %tmp(:,:,k) = nlfilter(rI(:,:,k),gridSize,fun);
        end
        tBlock = [cBlock,tBlock];
    end
    tBlock = single(tBlock);
    tBlock = squeeze(sum(tBlock,2));
end