function [b pth] = vTest3(image,D)
    fidx = find(image);
    pth = zeros(size(image));
    b = 0;
    if numel(fidx) > 2 & image(5) == 1 & sum(image(:)) <= 5
        %P = permn(fidx,numel(fidx));
        P = perms(fidx);
        for p = 1:size(P,1)
            %F = [];
            
            bFlag = 0;
            if P(p,1) ~= 5 & P(p,end) ~= 5;
                for c = 1:(size(P,2)-1)
                    bFlag = 1;
                    if D(P(p,c),P(p,c+1)) == 0
                        bFlag = 0;
                        break
                    end
                end
            end
            
            
            if bFlag
                b = 1;
                pth(P(p,:)) = 1:size(P,2);
                %P(p,:)
                %break
            end
            
            
        end
    end
    b;
end