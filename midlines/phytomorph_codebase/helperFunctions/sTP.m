function [R,IDX] = sTP(D,n,sz)

    v = 1:size(D,2);
    ci = nchoosek(v,n);
    R = zeros(size(D,1),size(ci,1));
    order = size(ci,2);
    dimy = prod(sz);
    MATSZ = dimy^order*ones(1,order);
    
    for e = 1:size(ci,1)
        for r = 1:size(ci,2)
            dump{r} = ci(e,r);
        end
        IDX(e) = sub2ind(MATSZ,dump{:});
    end
    
    for e = 1:size(D,1)
        tmp = ones(1,size(ci,1));
        for r = 1:size(ci,2)
            tmp = tmp .* D(e,ci(:,r));
        end
        R(e,:) = tmp;
    end
    
    
end

%{

[R,IDX] = sTP(rand(13,6),2,[2 3])

%}