function [T] = point2affine(P)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % construct basic affine transformation from point
    T = zeros(2,3,size(P,1));
    for pt = 1:size(P,1)        
        curPoint = [P(pt,2) P(pt,1)];
        T(:,:,pt) = [eye(2) curPoint'];
    end
end