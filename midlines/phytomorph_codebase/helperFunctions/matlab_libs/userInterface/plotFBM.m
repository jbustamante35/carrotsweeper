function [] = plotFBM(I,mapStack)
    imshow(I,[]);
    hold all;
    for i  = 1:size(mapStack,3)
        [x2 x1] = find(mapStack(:,:,i));
        plot(x1,x2,'.');
    end
end