function [Y] = unpackQbits(X,b)
    Y = sparse(size(X,1),size(X,2)*b);
    IDX = (1:size(X,1))';
    for byte = 1:size(X,2)
        for bit = 1:b
            tmp = find(bitget(X(:,byte),bit));
            ind = sub2ind(size(Y),IDX(tmp),(byte-1)*b + bit*ones(size(tmp)));
            Y(ind) = 1;
        end
    end
end