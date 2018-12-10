function [tipC] = getTip(rootCurve,SNIP)
    for i = 1:numel(rootCurve)
        % calculate curvature
        d1X1 = cwt(rootCurve{i}(:,1),SNIP,'gaus1');
        d1X2 = cwt(rootCurve{i}(:,2),SNIP,'gaus1');
        d2X1 = cwt(rootCurve{i}(:,1),SNIP,'gaus2');
        d2X2 = cwt(rootCurve{i}(:,2),SNIP,'gaus2');
        K = (d1X1.*d2X2 - d1X2.*d2X1).*(d1X1.^2 + d1X2.^2).^-3/2;
        K(1:2*SNIP) = 0;
        K(end-2*SNIP:end) = 0;
        [~,tipC(i)] = max(abs(K));
    end
end