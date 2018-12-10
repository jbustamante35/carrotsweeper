function [pointList,USE] = zoomGather_nonDB(fileList,BOX,toZoom)
    for e = 1:numel(fileList)
        I = imread(fileList{e});
        imshow(I,[]);
        title(num2str(e))
        USE{e} = questdlg('Use');
        xy = [];
        [sxy(:,1),sxy(:,2),V] = impixel(I);
        if toZoom
            for bx = 1:size(sxy,1)
                b = point2Box(sxy(bx,:),BOX);

                rpidx = b(1:2) <= 0;
                b(rpidx) = 1;

                rpidx = (b(1:2) + b(3:4)) > [size(I,2) size(I,1)];
                b(3:4) = b(3:4).*~rpidx + [size(I,2) size(I,1)].*rpidx;


                subI = imcrop(I,b);

                [tmp_xy(1),tmp_xy(2),V] = impixel(subI);

                sxy(bx,:) = tmp_xy + b(1:2);

            end
        end
        
        pointList{e} = sxy;
        
        close all
        imshow(I,[]);
        hold on
        plot(sxy(:,1),sxy(:,2),'g*');
        hold off
        drawnow
    end
    
end