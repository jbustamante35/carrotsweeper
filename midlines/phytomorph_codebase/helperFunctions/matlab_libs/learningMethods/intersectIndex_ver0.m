function [Mindex] = intersectIndex_ver0(Mindex,featureFiles)
    for e = 1:numel(featureFiles)
        load(featureFiles{e},'index');
        Mindex = intersect(Mindex,index);
    end
end