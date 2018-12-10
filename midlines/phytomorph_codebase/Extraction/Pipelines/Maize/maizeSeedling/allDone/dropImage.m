function [into,intoBlend,Z,wholeMask] = dropImage(into,boundingBox,toDrop,toDropMask,intoBlend)
    sz = size(into);
    
    UL = fliplr(round(boundingBox(1:2)));
    UL(UL <= 0) = 1;
    rec = [UL UL + [size(toDrop,1) size(toDrop,2)]-1];
    WID = round(fliplr(boundingBox(3:4))) + UL;
    
   
    wholeMask = zeros(size(into,1),size(into,2));
    
    if ~isempty(toDropMask)
        into(rec(1):rec(3),rec(2):rec(4),:) = toDrop.*toDropMask;
        wholeMask(rec(1):rec(3),rec(2):rec(4)) = toDropMask;
    else
        into(rec(1):rec(3),rec(2):rec(4),:) = toDrop;
    end
 
    
    
    intoBlend(rec(1):rec(3),rec(2):rec(4)) = toDropMask;
    
    
    Z = zeros(size(into,1),size(into,2));
    CP = [UL;WID];
    CP = round(mean(CP,1));
    Z(CP(1),CP(2)) = 1;
    
    
    Z = Z(1:sz(1),1:sz(2));
    
    into = into(1:sz(1),1:sz(2),:);
    intoBlend = intoBlend(1:sz(1),1:sz(2));
end