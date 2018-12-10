function [ktp] = projectTPToVforK(K,D,V,parameters,tp)
    for i = 1:size(tp,1)        
        ktp(:,i) = kernel1(D(:,1:size(tp,2)),tp(i,:),{parameters});
        %ktp(:,i) = center(ktp(:,i),K);
    end 
    ktp = center(ktp,K);
    ktp = (V'*ktp)';
end

function [tpc] = center(tp,K)
    n = size(K,1);
    Oh = ones(n,1);    
    O = Oh*Oh';
    KON = - K*Oh/n + O*K*Oh*n^-2;
    tpc = tp - O*tp*n^-1 + repmat(KON,[1 size(tp,2)]);
end


function [D] = kernel1(D,x,para)
    D = bsxfun(@minus,D,x);
    D = sum(D.*D,2);
    dist = D.^.5;
    D = D / (2 * para{1}^2);
    D = exp(-D);
    fidx = find(dist < para{1});
    D(fidx) = 0;
end