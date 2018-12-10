function [X,Y] = sampleSparkleSequence(FileList,sparkleSize,sparkleList,disp)
   
    % for each image
    for e = 1:numel(FileList)
        close all
        I = double(imread(FileList{e}))/255;
        I = imresize(I,sparkleSize);
        if e == 1
            X = zeros([size(I) numel(FileList)]);
        end
        X(:,:,:,e) = I;
        sz = size(sparkleList{e});
        Y(e,:) = reshape(sparkleList{e},[1 prod(sz)]);
        e
    end
end