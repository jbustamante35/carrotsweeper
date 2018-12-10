function [pt] = markFrame(T,BL,I)
    R = [];
    for c = 1:size(BL,3)
        R(:,:,c) = imfilter(double(T),double(BL(:,:,c)),'replicate');
    end
    [IDXP,IDX] = max(R,[],3);
    eT = imerode(T,strel('square',5));
    mIDX = IDX.*eT;
    %imshow(mIDX,[]);
    toSearchFor = [1 3 4 6 7 8 9 10];
    pt = [];
    for e = 1:numel(toSearchFor)
        ff = bwlarge(mIDX == toSearchFor(e));
        fidx = find(ff);
        [mm,midx] = max(IDXP(fidx));
        midx = find(IDXP(fidx) == mm);
        AL = [];
        [AL(:,1) AL(:,2)] = ind2sub(size(T),fidx(midx));
        pt(e,:) = mean(AL,1);
    end
    
    
    if ~isempty(I)
        imshow(I,[]);
        hold on
        plot(pt(:,2),pt(:,1),'ko')
        for pp = 1:size(pt,1)
            text(pt(pp,2),pt(pp,1),num2str(pp));
        end
        drawnow
        hold off
    end
end