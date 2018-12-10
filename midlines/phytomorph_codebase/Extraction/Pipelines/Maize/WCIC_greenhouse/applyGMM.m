function [MASK] = applyGMM(fileList,loadFunc,GMModel)
    % apply model
    for e = 1:numel(fileList)
        tmp = imread(fileList{e});
        tmp = loadFunc(double(tmp)/255);
        tmp = imfilter(tmp,fspecial('disk',11),'replicate');
        %tmp = decorrstretch(tmp);
        
        sz = size(tmp);
        tmpo = tmp;
        tmp = double(reshape(tmp,[size(tmp,1)*size(tmp,2) size(tmp,3)]));
        idx = cluster(GMModel,255*tmp);
        idx = reshape(idx,sz(1:2));
        MASK(:,:,e) = idx;
        %RGB = label2rgb(idx);
        %imshow([RGB tmpo],[]);
        %drawnow
    end
end