function straightMask = sampleStraighten(midline, carrotMask, vis)
try
    [DomainS, DomainG] = extendCarrotMidline(midline, [0 0], carrotMask);
    dsz = size(DomainG);
    
    straightMask = ba_interp2(double(carrotMask)/255,DomainS(:,2),DomainS(:,1));
    straightMask = reshape(straightMask,[dsz(1) dsz(2)]);
catch e
    fprintf(2, 'Error straightening mask\n%s\n', e.getReport);
    straightMask = [];
end

if vis
    imshow(straightMask, []);
    drawnow;
end

end