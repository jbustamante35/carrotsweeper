function [E] = mutualEntropy(X,Y)
    Z = [X,Y];

    [uqX] = generateUQbits(size(X,2));
    [uqY] = generateUQbits(size(Y,2));
    [uqZ] = generateUQbits(size(Z,2));
    O = ones(size(Z,1),1);
    for u = 1:size(uqX,1)
        %idxX = all(bsxfun(@eq,X,uqX(u,:)),2);
        idxX = O;
        for b = 1:size(X,2)
            idxX = idxX & (X(:,b) == uqX(u,b));
        end
        NX(u) = sum(idxX);
    end
    PX = NX*size(X,1)^-1;
    for u = 1:size(uqY,1)
        %idxY = all(bsxfun(@eq,Y,uqY(u,:)),2);
        idxY = O;
        for b = 1:size(Y,2)
            idxY = idxY & (Y(:,b) == uqY(u,b));
        end
        NY(u) = sum(idxY);
    end
    PY = NY*size(Y,1)^-1;
    
    for u = 1:size(uqZ,1)
        %idxZ = all(bsxfun(@eq,Z,uqZ(u,:)),2);
        idxZ = O;
        for b = 1:size(Z,2)
            idxZ = idxZ & (Z(:,b) == uqZ(u,b));
        end
        NZ(u) = sum(idxZ);
        idxX = uqZ(u,1:size(X,2));
        idxY = uqZ(u,size(X,2)+1:end);
        idxX = all(bsxfun(@eq,uqX,idxX),2);
        idxY = all(bsxfun(@eq,uqY,idxY),2);
        F(u) = PX(idxX)*PY(idxY);
    end
    PZ = NZ*size(Z,1)^-1;
    F = -log(PZ.*F.^-1);
    F(isinf(F)) = 0;
    E = -sum(F.*PZ);
end