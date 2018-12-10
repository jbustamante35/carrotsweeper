function [nearPT,farPT,positivePT] = generateSamplePoints_forBOX(pointMask,dilateValue,scaleValue,falseDistance,numberFalsePoints,numberPositivePoints)
    sz = size(pointMask);
    if dilateValue ~= 0
        if scaleValue ~= 1
            [px,py] = find(pointMask);
            pointMask = imresize(pointMask,scaleValue,'nearest');
            pointMask = zeros(size(pointMask));
            pointMask(round(py*scaleValue),round(px*scaleValue)) = 1;
            pointMask = pointMask > .5;
        end
        pointMask = imdilate(pointMask,strel('square',round(dilateValue*scaleValue)));
        if scaleValue ~= 1
            pointMask = imresize(pointMask,sz);
            pointMask = pointMask > .5;
        end
    end
   
   distT = bwdist(pointMask) < falseDistance & pointMask == 0;
   
   [nearPT(:,1),nearPT(:,2)] = find(distT);
   [posPT(:,1),posPT(:,2)] = find(pointMask);
   [farPT(:,1),farPT(:,2)] = find(distT==0 & pointMask == 0);
   
   
   r1 = randi(size(nearPT,1),[min(size(nearPT,1),numberFalsePoints(1)) 1]);
   r2 = randi(size(farPT,1),[min(size(farPT,1),numberFalsePoints(2)) 1]);
   
   if numberPositivePoints(1) > 0
        r3 = randi(size(posPT,1),[min(size(posPT,1),numberPositivePoints(1)) 1]);
   else
       r3 = 1:size(posPT,1);
   end
   
   nearPT = nearPT(r1,:);
   farPT = farPT(r2,:);
   %falsePT = [nearPT(r1,:);farPT(r2,:)];
   positivePT = [posPT(r3,:)];
   if any(farPT(:) < 0)
       here =1;
   end
end


%{
    [x,y,v] = impixel(I);
    pointMask = zeros(size(I,1),size(I,2));
    pointMask(y,x) = 1;
    [falsePT,positivePT] = generateSamplePoints_forBOX(pointMask,21,.25,300,[10 10],10);
    [falseBOX] = point2Box(falsePT,[100 100]);
    [positiveBOX] = point2Box(positivePT,[100 100]);
    sI = insertShape(I, 'Rectangle', falseBOX,'Color',{'red'},'LineWidth',3);
    sI = insertShape(sI, 'Rectangle', positiveBOX,'Color',{'green'},'LineWidth',3);
    imshow(sI,[]);
%}