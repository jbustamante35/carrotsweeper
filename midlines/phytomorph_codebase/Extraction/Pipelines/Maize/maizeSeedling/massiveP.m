function [LP] = massiveP(S)
    try
        tmp = logical(S > .8);
         % get the skeleton
        skeleton = bwmorph(tmp,'thin',inf);

        skel = skeleton;
        [basePoint] = getStemBasePoint(S,9,skel);
        [path,pathcost,midx,DP] = traceSeedlingSkeleton(skel,basePoint);
        LP = [DP(2,path{midx})',DP(1,path{midx})'];
    catch
        LP = [];
    end
end