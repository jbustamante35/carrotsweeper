function [R] = mask2cellBox(M,areaT,expand,disp,I,CL)
    redTape = M(:,:,1);
    redTape = bwareaopen(redTape,areaT);
    redTape = imclose(redTape,strel('line',400,90));
    redTape = imclose(redTape,strel('line',400,0));


    innerCell = imfill(redTape,'holes');
    innerCell = innerCell == 1 & redTape == 0;

    innerCell = bwareaopen(innerCell,areaT);

    R = regionprops(innerCell,'BoundingBox','Area');
    n = sum(count([R.Area]));

    innerCell = bwlarge(innerCell,n);
    R = regionprops(innerCell,'BoundingBox');
    for e = 1:numel(R)
        R(e).BoundingBox(1:2) = R(e).BoundingBox(1:2) - expand;
        R(e).BoundingBox(3:4) = R(e).BoundingBox(3:4) + 2*expand;
    end

    if disp
        out = I;
        for e = 1:size(M,3)
            out = flattenMaskOverlay(out,logical(M(:,:,e)),.3,CL{e});
        end

        imshow(out,[]);
        hold on 
        for b = 1:numel(R)
            rectangle('Position',R(b).BoundingBox,'EdgeColor','b','LineWidth',3);
        end
    end

end