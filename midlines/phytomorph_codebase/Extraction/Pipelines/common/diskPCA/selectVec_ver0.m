function [vec] = selectVec_ver0(vec,grp,e,k)
    vec = vec(:,grp(e,:)==k);
end