function [Z] = makeRingImage_org(oimg,img,grid)
    G1 = squeeze(grid(:,1,:));
    G2 = squeeze(grid(:,2,:));
    Z = zeros([size(oimg,2) size(oimg,1) size(img,2)]);
    for e1 = 1:size(G1,2)
        for e2 = 1:size(G2,1)
            for k = 1:size(img,2)
                try
                    x1p = round((G1(e2,e1)));
                    x2p = round((G2(e2,e1)));
                    Z(x2p,x1p,k) = img(e2,k,e1,1);
                catch

                end
            end
        end
    end    
    Z = (Z');
end