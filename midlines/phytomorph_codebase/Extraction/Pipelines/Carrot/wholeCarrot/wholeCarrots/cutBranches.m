function [img] = cutBranches(img,toRepair)

    pidx = find(imdilate(bwmorph(img,'branchpoints',1),ones(3)));
    img(pidx) = 0;
    
end