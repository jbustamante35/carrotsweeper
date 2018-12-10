function [F] = myT1(M)
    DIST = bwdist(M);
    DIST1 = bwdist(~M);
    F = DIST1 - DIST;
    for itr = 1:100
        F = F + .01*del2(M);
        %imshow(F,[]);
        %drawnow
        itr
    end
end