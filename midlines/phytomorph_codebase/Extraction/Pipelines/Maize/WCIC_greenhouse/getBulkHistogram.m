function [H] = getBulkHistogram(fileList)
    H = [];
    for e = 1:numel(fileList)
        I = imread(fileList{e});
        tmp = [];
        for k = 1:3
            tmp = [tmp imhist(I(:,:,k))'];
        end
        H = [H;tmp];
        fprintf(['Done with image:' num2str(e) '\n'])
    end
end