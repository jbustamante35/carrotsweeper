function [gamma] = liftCurve(curve,length,frame)
    for e = 1:size(curve.d,1)-length
        seg = curve(e:e+length);
        seg = bsxfun(@minus,seg,seg(1,:));
        
        tFrame = squeeze(frame(e));
        seg = (tFrame*seg')';
        gamma(e,:,:) = seg(:,1:2);
    end
end