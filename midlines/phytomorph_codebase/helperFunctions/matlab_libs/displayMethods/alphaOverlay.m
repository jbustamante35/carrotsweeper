function [h] = alphaOverlay(I,patchSet,h)
    if isempty(h)
        h = figure;
    end
    figure(h);
    imshow(I,[])
    for i = 1:numel(patchSet)
        hold on
        for j = 1:size(patchSet{i}.data,1)
            patch(patchSet{i}.data{j}(:,2),patchSet{i}.data{j}(:,1), patchSet{i}.color, 'FaceAlpha', 0.3)
        end
    end
end