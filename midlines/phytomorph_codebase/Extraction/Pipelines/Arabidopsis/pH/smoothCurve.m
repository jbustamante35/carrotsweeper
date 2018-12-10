function [y] = smoothCurve(dB)

    %f = spap2(7,3,1:size(dB{1},1),dB{1}');
    %y{1} = fnval(f,1:size(dB{1},1))';

    %y{1} = interp1(1:size(dB{1},1),dB{1},1:size(dB{1},1),'spline');
    f = fspecial('average',[31 1]);
    y{1} = [imfilter(dB{1}(:,1),f,'replicate') imfilter(dB{1}(:,2),f,'replicate')];
end