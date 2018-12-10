function [b] = vTest2(image)
    
    
    noHoles = (imfill(image,'holes') - image) == 1;
    noHoles = ~any(noHoles(:));
    
    Spur = (image - bwmorph(image,'spur',1)) == 1;
    Spur = any(Spur(:));
    
    b = noHoles & ~Spur;
end