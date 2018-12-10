function [vec,vecM] = sampleStraighten(midline,carrotMask,carrotImage)

    [DomainS,DomainG] = extendCarrotMidline(midline,[0 0],carrotMask);
    dsz = size(DomainG);
    vec = [];
    for k = 1:size(carrotImage,3)
        vec(:,k) = ba_interp2(double(carrotImage(:,:,k))/255,DomainS(:,2),DomainS(:,1));
    end

    vecM = ba_interp2(double(carrotMask)/255,DomainS(:,2),DomainS(:,1));


    vec = reshape(vec,[dsz(1) dsz(2) 3]);
    vecM = reshape(vecM,[dsz(1) dsz(2)]);
end