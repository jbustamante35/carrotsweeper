function [M] = connectPlant(M)
    % remove the main object - what is there is not a main object connected to the bottom
    mM = M - imclearborder(M);
    pM = M - mM;
    [d idx] = bwdist(mM);
    R = regionprops(logical(pM),'PixelIdxList');
    pM = mM;
    for e = 1:numel(R)
        dd = [];
        [d idx] = bwdist(pM);
        
        
        for p = 1:numel(R(e).PixelIdxList)
            dd(p) = d(R(e).PixelIdxList(p));
            map(p) = idx(R(e).PixelIdxList(p));
        end
        [J,midx] = min(dd);
        
        if J < 50
            [str(1) str(2)] = ind2sub(size(pM),map(midx));
            [stp(1) stp(2)] = ind2sub(size(pM),R(e).PixelIdxList(midx));
            DIS = dd(midx);
            X = round(linspace(str(1),stp(1),ceil(DIS)));
            Y = round(linspace(str(2),stp(2),ceil(DIS)));
            IDX = sub2ind(size(pM),X,Y);
            % connects via straight line
            pM(IDX) = 1;
            % adds the distant object
            pM(R(e).PixelIdxList) = 1;
        else
            
        end
        
        
    end
    M = pM;
end