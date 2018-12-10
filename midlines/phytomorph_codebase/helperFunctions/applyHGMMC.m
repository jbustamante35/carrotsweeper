function [idx] = applyHGMMC(I,hgmm,idx,currentLevel)

    if currentLevel == 1
        idx = cluster(I,);
    end


    for u = 1:numel(hgmm)
       
        for c = 1:hgmm(u).gmm.NumComponents

            fidx = find(idx==c);

            cluster(hgmm(u).gmm,I(fidx,:));
            
            subIDX = applyHGMMC(I(fidx,:),hgmm(u).nextLevel(c),idx(fidx));

            if u == 1
                nextLevel = zeros(size(idx,1),size(subIDX,2));
            end

            nextLevel(fidx,:) = subIDX;
        end
    end
    
end

function [newidx] = subCluster(I,gmm,idx)
    UQ = unique(idx);
    newidx = zeros(size(idx));
    for u = 1:numel(UQ)
        fidx = find(idx==UQ(u));
        newidx(fidx) = cluster(I,gmm);
    end
end