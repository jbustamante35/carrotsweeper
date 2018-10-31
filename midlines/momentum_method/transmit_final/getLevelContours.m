function [curve] = getLevelContours(I,levels)
    C = contourc(I,levels);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % get the curve structure
    str = 1;
    c = 1;
    clear curve
    while str < size(C,2)
        ed = str + C(2,str);
        curve(c).level = C(1,str);
        curve(c).data = C(:,str+1:ed);
        curve(c).length = size(curve(c).data,2);
        c = c + 1;
        str = ed + 1;
    end    
end