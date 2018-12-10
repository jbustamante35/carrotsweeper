function [mea,POS] = splitEnergyLevels(I,idx,distribution,MASK)

        %{
        for k = 1:3
            [d1(:,:,k) d2(:,:,k)] = gradient(I(:,:,k));
        end
        VEC = (d1(:,:,1).^2 + d1(:,:,2).^2 + d1(:,:,3).^2).^-.5;
        VEC = d1.*VEC;
        NOR = d2 - VEC.*dot(d2,VEC,3);
        NOR = (NOR(:,:,1).^2 + NOR(:,:,2).^2 + NOR(:,:,3).^2).^.5;
        TAN = (d1(:,:,1).^2 + d1(:,:,2).^2 + d1(:,:,3).^2).^.5;
        AREA = TAN.*NOR;

        mea = zeros(size(AREA));
        mea(idx==1) = AREA(idx == 1);
        mea = mea.*double(MASK);
        %}

        tmpI = I(:,:,1);
        [d1 d2] = gradient(tmpI);
        mea = (d1.^2 + d2.^2).^5;
        mea = mea.*double(MASK);
        
        fidx = find(MASK);
        vq = interp1(distribution(2,:),distribution(1,:),mea(fidx),'linear',.001);
        PROB = zeros(size(MASK));
        PROB(fidx) = vq;

        PROB = PROB / sum(PROB(:));


        [n1 n2] = ndgrid(1:size(PROB,1),1:size(PROB,2));


        POS = [PROB(:)'*n1(:) PROB(:)'*n2(:)];


end






















