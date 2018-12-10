function [I] = batchMatch(fileName,batchName)
    if ~isa(batchName,'double')
        for e = 1:numel(batchName)
            B(:,:,:,e) = imread(batchName{e});
        end
    else
        B = batchName;
    end
    if ischar(fileName)
        I = double(imread(fileName))/255;
    else
        I = fileName;
    end
    for i = 1:size(B,4)
        tmp(:,:,:,i) = imhistmatch(I,B(:,:,:,i),256);
    end
    tmp = mean(tmp,4); 
    I = mean(tmp,4);
end