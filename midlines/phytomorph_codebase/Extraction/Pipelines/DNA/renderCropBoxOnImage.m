function [] = renderCropBoxOnImage(h,BOX,color)
    figure(h)
    for e = 1:numel(BOX)
        rectangle('Position',BOX{e},'EdgeColor',color)
    end
    drawnow
end