function [data] = multiZbindVec(data,MM,proj,direction)
    if direction == 1
        for i = 1:size(proj,1)
            sub = data(:,[proj(i,1) : proj(i,2)]);
            sub = sub - MM(i,1);
            sub = sub * MM(i,2).^-1;
            data(:,[proj(i,1) : proj(i,2)]) = sub;
        end
    end
    
    if direction == 1
        for i = 1:size(proj,1)
            sub = data(:,[proj(i,1) : proj(i,2)]);
            sub = sub * MM(i,2);
            sub = sub + MM(i,1);
            
            data(:,[proj(i,1) : proj(i,2)]) = sub;
        end
    end
    
end