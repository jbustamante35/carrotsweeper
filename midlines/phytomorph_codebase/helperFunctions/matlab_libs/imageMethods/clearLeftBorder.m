function [B] = clearLeftBorder(B)
    B(1,:) = 0;
    B(end,:) = 0;
    B(:,end) = 0;

    CC = bwconncomp(B);
    B = zeros(size(B));
    for i = 1:CC.NumObjects
        [P1 P2] = ind2sub(size(B),CC.PixelIdxList{i});
        if any(P2==1)
            B(CC.PixelIdxList{i}) = 1;
        end
    end
end