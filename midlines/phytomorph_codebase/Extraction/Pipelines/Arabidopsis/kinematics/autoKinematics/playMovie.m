function [] = playMovie(imageStack)
    for e = 1:numel(imageStack)
        img = imread(imageStack{e});
        imshow(img,[]);
        drawnow
    end
end