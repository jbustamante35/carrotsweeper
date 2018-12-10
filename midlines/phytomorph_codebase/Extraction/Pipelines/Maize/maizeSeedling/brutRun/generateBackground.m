function [simBK] = generateBackground(vB1,hB1,vB2,hB2,sz)
    vB1 = repmat(vB1,[sz(1) 1 1]);
    vB2 = repmat(vB2,[1 sz(2) 1]);
end