function [vB1,hB1,vB2,hB2,TOT,BK] = generateBKsignals_wrapper(I,obj)

        sz = size(I);
        % reshape the data
        tmp1 = reshape(I,[prod(sz(1:2)) sz(3)]);

        [K,~,nl] = obj.cluster(tmp1);
        [K] = reshape(K,sz(1:2));
        [nl] = reshape(nl,[sz(1:2) size(nl,2)]);
        BK = K == 1;


        
        [vB1,hB1,vB2,hB2,TOT] = generateBKsignals(I,BK);
end