function [r] = lowMemCNN(block,net,sz,index,DISK)
    block = reshape(block,[sz 1 size(block,2)]);
    
     if nargin == 5
        block = ba_interp2(block,DISK(:,:,1),DISK(:,:,2));
    end
    
    [class prob] = classify(net,block);
    r = prob(:,index);
end