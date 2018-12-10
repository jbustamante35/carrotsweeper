function [IDX BOX] = generateImageIDXPatches(SZ,B,N)
    R = (B+1):(SZ(1)-B);
    C = (B+1):(SZ(2)-B);
    LENR = floor(numel(R)/N(1));
    LENC = floor(numel(C)/N(2));
    REMR = rem(numel(R),N(1));
    REMC = rem(numel(C),N(2));
    cnt = 1;
    str1 = 1;
    for n1 = 1:N(1)
        stp1 = str1 + LENR - 1;
        
        if n1 == N(1)
            stp1 = stp1 + REMR;
        end
        str2 = 1;
        for n2 = 1:N(2)
            stp2 = str2 + LENC - 1;
            if n2 == N(2)
                stp2 = stp2 + REMC;
            end
            MASK = zeros(SZ);
            MASK(R(str1):R(stp1),C(str2):C(stp2)) = 1;
            RE = regionprops(MASK==1,'BoundingBox');
            BOX(cnt).BoundingBox = RE(1).BoundingBox;
            IDX(cnt).R = [R(str1):R(stp1)];
            IDX(cnt).C = [C(str2):C(stp2)];
            cnt = cnt + 1;
            str2 = stp2 + 1;
        end
        str1 = stp1 + 1;
    end
end