function [I DIFF] = fixImage(I,brightI,colorI,n)
    oI = I;
    for e = 1:n
        HSV = rgb2hsv(I);
        HSV(:,:,3) = imhistmatch(HSV(:,:,3),brightI,10);
        I = hsv2rgb(HSV);
        I = imhistmatch(I,colorI,255);
        I = mean(cat(4,I,oI),4);
        DIFF(e) = norm(I(:) - oI(:));
    end
end