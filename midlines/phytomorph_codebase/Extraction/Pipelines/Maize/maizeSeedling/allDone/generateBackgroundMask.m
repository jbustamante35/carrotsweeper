function [MSK] = generateBackgroundMask(I,TOT)
    diffI = I-TOT;
    diffI = sum(diffI.*diffI,3).^.5;
    bx = linspace(0,1,300);
    v = hist(diffI(:),bx);
    v = v / sum(v);
    plot(bx,v);

    TH = graythresh(diffI);
    MSK = diffI > TH*1.05;

    MSK = imclose(MSK,strel('disk',5,0));

    top = MSK(1,:);
    bot = MSK(end,:);
    MSK(1,:) = 1;
    MSK(end,:) = 1;
    MSK = imfill(MSK,'holes');
    MSK(1,:) = top;
    MSK(end,:) = bot;
end