function [d p] = whiten(d,p)
    
    if nargin == 1
        p = std(d,1,1);
    end
    
    for i = 1:size(d,1)
        d(i,:) = d(i,:).*p.^-1;
    end
end