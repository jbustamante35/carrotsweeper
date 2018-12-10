function [u edge] = measureUCG(pt)
    
    p = nchoosek(1:size(pt,1),2);
    for e = 1:size(p,1)
        edge(:,:,e) = pt(p(e,:),:);
        u(e) = norm(pt(p(e,1),:)-pt(p(e,2),:));
    end
end