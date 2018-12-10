function [stack] = grabStack(FileList,str,skp,stp,BOX)
    vec = str:skp:stp;
    cnt = 1;
    if nargin == 4
        stack= imread(FileList{1});
    else
        BOX = round(BOX);
        IDX = {[BOX(2) 1 BOX(2) + BOX(4)] [BOX(1) 1 BOX(1)+BOX(3)]};
        stack(:,:,:,cnt) = imread(FileList{1},'PixelRegion',IDX);
    end
    stack = zeros([size(stack) numel(vec)]);
    for e = 1:numel(vec)
        if nargin == 4
            stack(:,:,:,cnt) = imread(FileList{e});
        else
            IDX = {[BOX(2) 1 BOX(2) + BOX(4)] [BOX(1) 1 BOX(1)+BOX(3)]};
            stack(:,:,:,cnt) = imread(FileList{e},'PixelRegion',IDX);
        end
        cnt = cnt + 1;
        cnt
    end
end