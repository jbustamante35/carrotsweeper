function [sig var] = subMBL(sig,fFrames)
    for e = 1:size(sig,2)
        sig(:,e) = subBL(sig(:,e),fFrames);
    end
end