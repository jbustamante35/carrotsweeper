function [patch] = map2patch(mapStack,direction,level)
    for i = 1:size(mapStack,3)
        if isempty(level)
            level = graythresh(mapStack(:,:,i));
        end
        switch direction
            case 'l'
                b = mapStack(:,:,i) < level;
            case 'g'
                b = mapStack(:,:,i) > level;
        end
        patch = bwboundaries(b);
    end
end