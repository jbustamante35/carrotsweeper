function [stack] = pushPop(FileList,ptr,stack,N,BOX)
    for e = 1:N
        stack(:,:,:,1) = [];
        if nargin == 4
            tmp = imread(FileList{1});
        else
            BOX = round(BOX);
             IDX = {[BOX(2) 1 BOX(2) + BOX(4)] [BOX(1) 1 BOX(1)+BOX(3)]};
            tmp = double(imread(FileList{ptr},'PixelRegion',IDX));
            %tmp = rgb2gray(tmp);
        end
        stack = cat(4,stack,tmp);
    end
end