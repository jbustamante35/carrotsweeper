function [b] = vTest(image,P,testEdge,holeTest)
    
    s = [image(1:end-1,1);image(end,1:end-1)';image(2:end,end);image(1,2:end)';];
    
    
    noHoles = (imfill(image,'holes') - image) == 1;
    noHoles = ~any(noHoles(:));
    
    Spur = (image - bwmorph(image,'spur',1)) == 1;
    Spur = any(Spur(:));
    
    R = bwconncomp(logical(image));
    
    image = im2colF(image,[3 3],[1 1]);
    fidx = find(image(5,:)==1);
    image = (2*ones(1,9)).^fliplr((0:8))*image(:,fidx)+1;
    
    if nargin ~= 2
        noHoles = 1;
    end
    if nargin == 2
        b = all(P(image)) & sum(s) == 2 & R.NumObjects == 1 & noHoles & ~Spur;
    else
        b = all(P(image)) & R.NumObjects == 1 & noHoles & Spur;
    end
end