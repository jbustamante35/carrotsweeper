function [] = displayModel(affine,comp,dataID)
    for s = 1:numel(affine)
        for t = 1:size(affine{s},3)
            for e = 1:size(comp{s},3)
                temp{s}(:,t,e) = affine{s}(:,:,t)*[comp{s}(:,t,e);1];
            end
        end
    end
    CL = {'r' 'b' 'g'};
    for t = 1:size(temp{1},2)
        for s = 1:numel(temp)
            midline = squeeze(temp{s}(:,t,:));
            s
            plot(midline(1,:),midline(2,:),CL{s});
            hold all
            quiver(affine{s}(1,3,t),affine{s}(2,3,t),affine{s}(1,1,t),affine{s}(2,1,t),100);
            axis equal
        end
        drawnow
        hold off
    end
end