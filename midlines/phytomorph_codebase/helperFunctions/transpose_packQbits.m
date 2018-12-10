function [Y] = transpose_packQbits(X,b)
    cnt = 1;
    for col = 1:size(X,2)
        tmp = zeros(1,size(X,1));
        for bit = 1:b
            tmp = bitset(tmp,b,bitget(X(:,col),b)');
        end
        Y(col,:) = squeezeQbits(tmp,8,3,'uint8');
        col
    end
end