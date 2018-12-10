function [] = renderBoxesOnImage(I,box,color,width)
    if ischar(I)
        I = double(imread(I))/255;
    end

    imshow(I,[]);
    drawnow
    for b = 1:numel(box)
        rectangle('Position',box{b},'EdgeColor',color,'LineWidth',width);
    end
end