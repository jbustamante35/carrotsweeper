function [box] = orderRedCropBoxes(box)
    try
        pointList = [];
        idList = (1:9)';

        for c = 1:numel(box)
            pointList(c,:) = box{c}(1:2);
        end

        nidx = [];
        for row = 1:3
            distance = sum(pointList.*pointList,2);
            % find the upper left corner
            [~,midx] = min(distance);
            % find the first row
            rowDistance = pointList(:,2) - pointList(midx,2);
            [~,sidx] = sort(rowDistance);

            [~,sidx1] = sort(pointList(sidx(1:3),1));


            nidx = [nidx;idList(sidx(sidx1(1:3)))];


            pointList(sidx(1:3),:) = [];
            idList(sidx(1:3)) = [];

        end
        
        
        box = box(nidx);
        %{
        for e = 1:size(pointList,1)
            plot(pointList(nidx(e),1),pointList(nidx(e),2),'r*')
            text(pointList(nidx(e),1),pointList(nidx(e),2),num2str(e))
            hold on
            waitforbuttonpress
        end
        %}
        
        
    catch
    end
end