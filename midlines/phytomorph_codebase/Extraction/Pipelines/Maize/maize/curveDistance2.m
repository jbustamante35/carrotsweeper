function [kD] = curveDistance2(target_fn,source_fn,lengthTarget,lengthSource,sampleReg)
    

    Ltarget = [0 cumsum(ones(1,size(sampleReg,2)))];
    Ltarget = Ltarget.*Ltarget(end)^-1;
    Ltarget = Ltarget*lengthTarget;
    
    Lsource = [zeros(size(sampleReg,1),1) cumsum(sampleReg,2)];
    Lsource = bsxfun(@times,Lsource,Lsource(:,end).^-1);
    Lsource = Lsource*lengthSource;
    
    
    % calculate curvature
    d_target = fnder(target_fn,1);
    d_source = fnder(source_fn,1);
    
    d_source = fnval(d_source,Lsource);
    d_target = fnval(d_target,Ltarget);
  
    D_source = sum(d_source.*d_source,1).^-.5;
    D_target = sum(d_target.*d_target,1).^-.5;
    
    d_source = bsxfun(@times,d_source,D_source);
    d_target = bsxfun(@times,d_target,D_target);
    d_target = [-d_target(2,:);d_target(1,:)];
    d_target = reshape(d_target,[size(d_target,1) 1 size(d_target,2)]);
    kD = bsxfun(@times,d_source,d_target);
    kD = squeeze(sum(sum(kD,1).^.5,3));
end