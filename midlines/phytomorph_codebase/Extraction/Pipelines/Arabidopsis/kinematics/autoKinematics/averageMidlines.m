function [AmidlineM] = averageMidlines(midlineM)
    for m = 1:numel(midlineM)
        path = midlineM{m};
        path = imfilter(path,fspecial('average',[100 1]),'replicate');
        path = arcLength(path,'arcLen');
        NUM(m) = size(path,1);
    end
    np = min(NUM);
    midlineStack = [];
    for m = 1:numel(midlineM)
        path = midlineM{m};
        path = imfilter(path,fspecial('average',[100 1]),'replicate');
        path = arcLength(path,'arcLen');
        midlineStack = cat(3,midlineStack,path(1:np,:));
    end
    AmidlineM = mean(midlineStack,3);
end