function [K_n,V,tp] = KPCA(X,N,sigma)
    K_n = K_nc(X,{sigma},@(x,y,p)kernel1(x,y,p));
    K_c = centerK(K_n);
    [V,D] = eigs(K_c,N);
    %{
    for e = 1:size(V,2)
        A = (size(X,1)*D(e,e))^-.5;
        V(:,e) = V(:,e)*A;%^.5; 
    end
    %}
    tp = (V'*K_c)';
end

function [Kc] = centerK(K)
    n = size(K,1);
    O = ones(n,1)*ones(1,n)/n;
    Kc = K - O*K - K*O + O*K*O;
end

function [K] = K_nc(D,parameters,kFunc)
    K = zeros(size(D,1));
    for i = 1:size(D,1)
        K(i,:) = kFunc(D,D(i,:),parameters);
    end
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