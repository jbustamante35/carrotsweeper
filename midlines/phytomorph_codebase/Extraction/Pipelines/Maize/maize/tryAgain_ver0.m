for e = 1:numel(SET)
    % read the image
    I = double(imread(SET{e}{1}));
    % complement the image
    I = imcomplement(I);
    % erode the image
    eI = imerode(I,strel('disk',51));
    % reconstruct
    BK = imreconstruct(eI,I);
    % subtract off background
    I = I - BK;
    imshow(I,[]);
    
    
end