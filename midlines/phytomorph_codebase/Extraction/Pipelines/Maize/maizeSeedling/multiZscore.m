function [data] = multiZscore(data,mu,sigma,proj,direction)
    if direction == 1
        for i = 1:size(proj,1)
            sub = data(:,[proj(i,1) : proj(i,2)]);
            sub = bsxfun(@minus,sub,mu(i,:));
            sub = bsxfun(@times,sub,sigma(i,:).^-1);
            data(:,[proj(i,1) : proj(i,2)]) = sub;
        end
    end
    
    if direction == 1
        for i = 1:size(proj,1)
            sub = data(:,[proj(i,1) : proj(i,2)]);
            sub = bsxfun(@times,sub,sigma(i,:));
            sub = bsxfun(@plus,sub,mu(i,:));
            data(:,[proj(i,1) : proj(i,2)]) = sub;
        end
    end
    
end