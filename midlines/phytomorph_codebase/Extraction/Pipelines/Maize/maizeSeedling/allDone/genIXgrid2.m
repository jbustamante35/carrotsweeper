function [P,pSZ] = genIXgrid2(imageSZ,skipValue,padValue,cleanBorder)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create sample point(s) grid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L1 = (imageSZ(1)-padValue(1)) - (padValue(1)+1);
    L2 = (imageSZ(2)-padValue(2)) - (padValue(2)+1);
    P1 = round(L1 * skipValue(1)^-1);
    P2 = round(L2 * skipValue(2)^-1);
    
    
    [d1,d2] = ndgrid(linspace((padValue(1)+1),(imageSZ(1)-padValue(1)),P1),...
                     linspace((padValue(2)+1),(imageSZ(2)-padValue(2)),P2));
                 
                 
    pSZ = size(d1);
    P = [d1(:),d2(:)];
    if nargin == 4
        rm1 = (P(:,1) < cleanBorder(1) | P(:,1) >= (imageSZ(1) - cleanBorder(1)));
        rm2 = (P(:,2) < cleanBorder(2) | P(:,2) >= (imageSZ(2) - cleanBorder(2)));
        P(rm1|rm2,:) = [];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end