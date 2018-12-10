function [Z] = makeRingImage(img)

    Rmin = 20;
    EXTRA = 0;
    Rmax = Rmin + size(img,3) + EXTRA;
    [g1 g2] = ndgrid(linspace(Rmax,Rmin,size(img,3)),linspace(-pi,pi,size(img,1)));
    G1 = g1.*cos(g2-pi/2);
    G2 = g1.*sin(g2-pi/2);
    
    SZ = 3*Rmax;
    Z = zeros([SZ SZ size(img,2)]);
    for e1 = 1:size(G1,2)
        for e2 = 1:size(G2,1)
            for k = 1:size(img,2)
                x1p = round(G1(e2,e1) + SZ/2) + 1;
                x2p = round(G2(e2,e1) + SZ/2) + 1;
                Z(x2p,x1p,k) = img(e1,k,e2,1);
            end
        end
    end
end