function  [tmpMaskBLUE] = generateBlueBoxMarks(maskSZ,pointListM,scale,da)
 try
        for pl = 1:size(pointListM,1)
            pointList = squeeze(pointListM(pl,:,:));
            pointList = fliplr(pointList');
            
            PUSH = 7;
            if pl == 1
                pointList = bsxfun(@plus,pointList,[-PUSH PUSH]);
            elseif pl == 2
                pointList = bsxfun(@plus,pointList,[PUSH -PUSH]);
            elseif pl == 3
                pointList = bsxfun(@plus,pointList,[-PUSH -PUSH]);
            elseif pl == 4
                pointList = bsxfun(@plus,pointList,[PUSH PUSH]);
            end
            
            pointList = round(bsxfun(@times,pointList,scale));
            msk = zeros([maskSZ]);
            for e = 1:size(pointList,1)
                msk(pointList(e,1),pointList(e,2)) = 1;
            end
            msk = imdilate(msk,strel('disk',da,0));
            tmpMaskBLUE(:,:,pl) = msk;
        end
    catch ME
        ME
    end
end