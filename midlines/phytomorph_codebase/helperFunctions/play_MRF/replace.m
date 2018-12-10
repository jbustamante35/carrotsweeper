function [of1 of2 d1 d2 new_image] = replace(image,in_filter,baseFilter,f1,f2,T,P,target,selectPoint)
    zi = image;
    zi;
    zi(5) = 0;
    [p1 p2] = find(zi);
    p1 = p1 - 2;
    p2 = p2 - 2;
    d1 = [p1(1) p2(1)];
    d2 = [p1(2) p2(2)];
    
    image = im2colF(image,[3 3],[1 1]);
    idx = (2*ones(1,9)).^fliplr((0:8))*image+1;
    %idx = find(find(P) == idx); % OLD INDEX FOR TIGHT CASE
    new_pot = baseFilter(:,idx).*in_filter;
    fidx = find(new_pot);
    %new_image_number = fidx(1);
    % sidx = randi(numel(fidx),1); % FOR RANDOM GENERATE SHAPE
    
    
    for s = 1:numel(fidx)
        value(s) = T(fidx(s),:)*target(selectPoint);
    end
    [J sidx] = max(value);
    
    new_image_number = fidx(sidx);
    
    
    idx;
    
    new_image_number;
    
    new_image = T(new_image_number,:);
    of1 = f1(:,idx,new_image_number);
    of2 = f2(:,idx,new_image_number);
end