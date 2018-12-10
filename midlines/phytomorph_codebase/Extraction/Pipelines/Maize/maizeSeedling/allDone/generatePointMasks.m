function [msk] = generatePointMasks(maskSZ,pointList,scale,da)
    try
        pointList = round(bsxfun(@times,pointList,scale));
        msk = zeros([maskSZ size(pointList,1)]);
        for e = 1:size(pointList,1)
            msk(pointList(e,1),pointList(e,2),e) = 1;
            msk(:,:,e) = imdilate(msk(:,:,e),strel('disk',da,0));
        end
    catch ME
        ME
    end
end