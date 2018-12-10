function [P,pSZ] = genIXgrid(imageSZ,skipValue,padValue,cleanBorder)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create sample point(s) grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [d1,d2] = ndgrid((padValue(1)+1):skipValue(1):(imageSZ(1)-padValue(1)),...
                     (padValue(2)+1):skipValue(2):(imageSZ(2)-padValue(2)));
    pSZ = size(d1);
    P = [d1(:),d2(:)];
    if nargin == 4
        rm1 = (P(:,1) < cleanBorder(1) | P(:,1) >= (imageSZ(1) - cleanBorder(1)));
        rm2 = (P(:,2) < cleanBorder(2) | P(:,2) >= (imageSZ(2) - cleanBorder(2)));
        P(rm1|rm2,:) = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end