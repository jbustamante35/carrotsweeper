function [p] = funcG(pts)
    p = mean(diff(pts,1,2),2);
    %p = p / norm(p);
end