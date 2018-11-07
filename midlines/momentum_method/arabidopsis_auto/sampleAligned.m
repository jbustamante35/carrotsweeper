function [X] = sampleAligned(snip1,spaceFraction,timeFraction)
    [n1 n2] = ndgrid(linspace(1,spaceFraction*size(snip1,2),size(snip1,2)),...
            linspace(1,timeFraction*size(snip1,3),size(snip1,3)));
    X1 = interp2(squeeze(snip1(1,:,:)),n1,n2);
    X2 = interp2(squeeze(snip1(2,:,:)),n1,n2);
    X = cat(3,X1,X2);
    X = permute(X,[3 1 2]);
end