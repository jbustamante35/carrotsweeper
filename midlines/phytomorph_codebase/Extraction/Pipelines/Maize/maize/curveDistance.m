function [kD] = curveDistance(fn,totalLength,sampleReg,baseLine)
    L = [zeros(size(sampleReg,1),1) cumsum(sampleReg,2)];
    L = bsxfun(@times,L,L(:,end).^-1);
    L = L*totalLength;
    % calculate curvature
    d1 = fnder(fn,1);
    d2 = fnder(fn,2);
    d1 = fnval(d1,L);
    d2 = fnval(d2,L);
    K = squeeze((d1(1,:,:).*d2(2,:,:) - d1(2,:,:).*d2(1,:,:)).*(d1(1,:,:).^2 + d1(2,:,:).^2).^-3/2);
    K = cumsum(K,2);
    kD = bsxfun(@minus,K,baseLine);
    kD = sum(kD.*kD,2).^.5;
end