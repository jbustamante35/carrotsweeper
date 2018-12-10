function [dP] = measureDisplacement(contour,tipIndex)
    for e = 1:numel(contour)
        path(e,:) = contour{e}(tipIndex(e),:);
    end
    dP = diff(path,1,1);
    dP = sum(dP.*dP,2).^.5;
    dP = [0;cumsum(dP)];
end