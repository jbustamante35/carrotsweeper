function [] = generateAndPushThumbNail(I)
    c = size(I,3);
    I = imresize(I,[50 50 c]);
    
end