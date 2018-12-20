function [R,carrotImage,carrotMask,midline,contour,straightRGB,straightMSK] = cropTraceStraightenSample(I,CARROT_SQUARE,ORIN,qrMSG)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % trace midlines and contours
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp = 1;
    if disp
        image(I);
        axis equal
        hold on
        axis off
    end
    midline = {};
    contour = {};
    carrotMask = {};
    carrotImage = {};
    straightRGB = {};
    straightMSK = {};
    for b = 1:numel(CARROT_SQUARE)
        if ~isempty(qrMSG{b})
            try
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make box straight wrt vertical camera line
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [newV] = min(CARROT_SQUARE{b}(2:3,1));
                CARROT_SQUARE{b}(2:3,1) = newV;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get crop box
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                MSK = (poly2mask(CARROT_SQUARE{b}(:,1),CARROT_SQUARE{b}(:,2),size(I,1),size(I,2)));
                R{b} = regionprops(MSK,'BoundingBox');
                carrotImage{b} = imcrop(I,R{b}(1).BoundingBox);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if flip needed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~ORIN(b)
                    carrotImage{b} = flip(carrotImage{b},2);
                end



                HSV = rgb2hsv(carrotImage{b});
                LAB = rgb2lab(carrotImage{b});
                carrotMask{b} = (LAB(:,:,1) < 65);
                carrotMask{b} = bwlarge(carrotMask{b});
                carrotMask{b} = imfill(carrotMask{b},'holes');
                sig = sum(carrotMask{b},1);
                fidx = find(sig > 30);
                
                R{b}.BoundingBox(3) = R{b}.BoundingBox(3) + fidx(1) - 1;
                R{b}.BoundingBox(3) = R{b}.BoundingBox(3) - fidx(1) - 1;
                
                
                % recrop
                carrotImage{b} = imcrop(I,R{b}(1).BoundingBox);
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if flip needed
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~ORIN(b)
                    carrotImage{b} = flip(carrotImage{b},2);
                end
                
                
                
                HSV = rgb2hsv(carrotImage{b});
                  LAB = rgb2lab(carrotImage{b});
                carrotMask{b} = ~(HSV(:,:,1) > .1);
                carrotMask{b} = (LAB(:,:,1) < 65);
                carrotMask{b} = imfill(carrotMask{b},'holes');
                
                
                saveVec = carrotMask{b}(:,1);
                carrotMask{b}(:,1) = 0;
                carrotMask{b} = imclearborder(carrotMask{b});
                carrotMask{b} = bwlarge(carrotMask{b});
                saveVec(1:100) = 0;
                saveVec(end-100:end) = 0;
                saveVec = bwlarge(saveVec);
                saveVec = carrotMask{b}(:,2);
                carrotMask{b}(:,1) = saveVec;
                cidx = find(carrotMask{b}(:,1));
                carrotMask{b}(cidx(1):cidx(end),1) = 1;
                carrotMask{b} = imfill(carrotMask{b},'holes');
                carrotMask{b} = bwlarge(carrotMask{b});
               
                [midline{b},contour{b}] = getMidlineAndContour(carrotMask{b},ORIN(b));


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % sample straight midline
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [straightRGB{b},straightMSK{b}] = sampleStraighten(midline{b},flip(carrotMask{b},2),flip(carrotImage{b},2));



                midlineDISPLAY = bsxfun(@plus,midline{b},R{b}(1).BoundingBox(1:2));
                contourDISPLAY = bsxfun(@plus,contour{b},R{b}(1).BoundingBox(1:2));


                if disp
                    plot(midlineDISPLAY(:,1),midlineDISPLAY(:,2),'g')
                    plot(contourDISPLAY(:,1),contourDISPLAY(:,2),'m')
                    drawnow
                end

            catch ME

                getReport(ME)
            end
        end
    end
end