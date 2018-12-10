function [Y] = randBitSample(X,bitsToSample,bitMAX)
    Y = zeros(size(X,1),1,'uint8');
    RND = randi(size(X,2),bitsToSample,1);
    bitCNT = 1;
    for e = 1:numel(RND)
        b = randi(bitMAX,1,1);
        Y = bitset(Y,bitCNT,bitget(X(:,RND(e)),b));
        bitCNT = bitCNT + 1;
    end
end