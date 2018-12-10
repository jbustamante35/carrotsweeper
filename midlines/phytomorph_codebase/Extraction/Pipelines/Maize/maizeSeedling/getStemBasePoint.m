function [basePoint] = getStemBasePoint(MASK,SNIP,skeleton)
    % init basePoint to empty
    basePoint = [];
    % SNIP off some stem
    baseMASK = sum(MASK((end-SNIP):end,:),1);
    fidx = find(baseMASK==max(baseMASK));
    if ~isempty(fidx)
        % mean along the stem snip
        basePoint(1) = mean(fidx);
        % set the basepoint 1 to the size of the mask
        basePoint(2) = size(MASK,1);
        % find skeleton
        [x,y] = find(skeleton);
        % stack the skeleton points for snap to
        DP = [x y]';
        % snap the base point to the skeleton
        [idx(1)] = snapTo(DP',[fliplr(basePoint)]);
        basePoint = DP(:,idx(1))';
    end
end