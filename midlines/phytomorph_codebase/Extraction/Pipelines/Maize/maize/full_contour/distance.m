function [D] = distance(source,target)
    rot = source.orientation(:,1)'*target.orientation'
    bend = mean(sum(source.data - target.data).^2,1).^.5);
    imgD = sum(abs(source.sam(:) - target.sam(:)));
    D(1) = rot;
    D(2) = bend;
    D(2) = imgD;
end