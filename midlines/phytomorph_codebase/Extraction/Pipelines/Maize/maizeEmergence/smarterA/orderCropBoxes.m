function [boundingBox,centerPoints,LABELS] = orderCropBoxes(MASK,boundingBox,centerPoints,I,disp)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % map centers to cells
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(['Starting circle mapping process \n']);tm = clock;
    [sortedCenters,LABELS] = getMap_ver1(logical(MASK));
    if disp
       imshow(I,[]);
       hold on
    end
    for e = 1:size(sortedCenters,1)
        delta = bsxfun(@minus,centerPoints,sortedCenters(e,:));
        [J,sidx(e)] = min(sum(delta.*delta,2));
    end
    boundingBox = boundingBox(sidx);
    centerPoints = centerPoints(sidx,:);
    if disp
        for e = 1:size(centerPoints,1)
            rectangle('Position',boundingBox{e},'EdgeColor','y')
            text(centerPoints(e,1)+200,centerPoints(e,2)+200,num2str(e),'Background','w');
            hold off
            drawnow
        end
    end
    fprintf(['Ending circle mapping process:' num2str(etime(clock,tm)) '\n']);
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    
end