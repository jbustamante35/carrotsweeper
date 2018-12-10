function [V] = vSm(V,idx,toState)
    switch toState
        case 'mat'
            M = zeros(size(idx.toMat,1).^.5);
            M(idx.toMat(:,1)) = V(idx.toMat(:,2));
            V = M;            
        case 'vec'
            V = V(idx.toVec)';
    end
end