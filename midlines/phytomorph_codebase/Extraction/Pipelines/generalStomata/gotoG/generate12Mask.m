function [myMask1,myMask2,totalMask] = generate12Mask(maskTree,thresholdValues,sMSK)


        FL = [[maskTree.Area]' [maskTree.Eccentricity]' [maskTree.MeanIntensity]' [maskTree.DistanceToEdge]'];
        for e = 1:size(FL,1)
            if (FL(e,4) > thresholdValues(4))
                threshValue1(e) = thresholdValues(1);
                threshValue2(e) = thresholdValues(2);
                threshValue3(e) = thresholdValues(3);
            else
                threshValue1(e) = thresholdValues(5);
                threshValue2(e) = thresholdValues(6);
                threshValue3(e) = thresholdValues(7);
            end
        end
        
        filter = find((FL(:,1) > threshValue1') & ...
                      (FL(:,2) > threshValue2') & ...
                      (FL(:,3) > threshValue3'));
                  
                  
        myMask1 = zeros(size(sMSK));
        for o = 1:numel(filter)
            try
                myMask1(maskTree(filter(o)).PixelIdxList) = 1;
                myMask1 = imfill(myMask1,'holes');
            catch ME
                ME;
            end
        end



        RS = regionprops(logical(myMask1),'MajorAxisLength');
        singleLength = mean([RS.MajorAxisLength]);



        myMask2 = zeros(size(myMask1));
        filter = find(round([maskTree.MajorAxisLength]' * singleLength^-1)==2 & (FL(:,3) > .6));
        for o = 1:numel(filter)
            myMask2(maskTree(filter(o)).PixelIdxList) = 1;
        end
        myMask2 = zeros(size(myMask1));
        for r = 1:4
            orginalBorder{r} = myMask1(:,1);
            myMask1(:,1) = 1;
            myMask1 = imrotate(myMask1,90);
        end
        holeFilled = imfill(myMask1,'holes');
        holeFilled = myMask1 == 0 & holeFilled == 1;
        holeFilledBK = bwlarge(holeFilled);
        holeFilled = holeFilledBK == 0 & holeFilled == 1;
        myMask1 = myMask1 + holeFilled;
        for r = 1:4
            myMask1(:,1) = orginalBorder{r};
            myMask1 = imrotate(myMask1,90);
        end

        for r = 1:4
            orginalBorder{r} = myMask2(:,1);
            myMask2(:,1) = 1;
            myMask2 = imrotate(myMask2,90);
        end
        holeFilled = imfill(myMask2,'holes');
        holeFilled = myMask2 == 0 & holeFilled == 1;
        holeFilledBK = bwlarge(holeFilled);
        holeFilled = holeFilledBK == 0 & holeFilled == 1;
        myMask2 = myMask2 + holeFilled;
        for r = 1:4
            myMask2(:,1) = orginalBorder{r};
            myMask2 = imrotate(myMask2,90);
        end

        % upgrade single to doubles
        RS1 = regionprops(logical(myMask1),'PixelIdxList');
        RS2 = regionprops(logical(myMask2),'PixelIdxList');
        myMask1 = zeros(size(myMask1));
        myMask2 = zeros(size(myMask2));

        toRM = [];
        for r = 1:numel(RS1)
            for s = 1:numel(RS2)
                if ~isempty(intersect(RS2(s).PixelIdxList,RS1(r).PixelIdxList))
                    toRM = [toRM r];
                end
            end
        end
        RS1(toRM) = [];
        for r = 1:numel(RS1)
            myMask1(RS1(r).PixelIdxList) = 1;
        end
        for r = 1:numel(RS2)
            myMask2(RS2(r).PixelIdxList) = 1;
        end

        myMask1 = bwareaopen(myMask1,40);
        myMask2 = bwareaopen(myMask2,40);

        Rclick = regionprops(logical(sMSK),'PixelIdxList');
        totalMask = logical(myMask1 + myMask2);
        Rtotal = regionprops(logical(totalMask),'PixelIdxList');

        sMSK = logical(sMSK);
        totalMask = logical(totalMask);

end