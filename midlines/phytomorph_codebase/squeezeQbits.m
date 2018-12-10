function [Y] = squeezeQbits(X,bitSequence,LOG_BYTE_SIZE,TYPE)
    
    if numel(bitSequence) == 1
        bitSequence = bitSequence*ones(1,size(X,2));
    end
    
    NUMBER_BYTES = ceil(sum(bitSequence)/(2^LOG_BYTE_SIZE));
    Y = zeros(size(X,1),NUMBER_BYTES,TYPE);
    
    newColumn = 1;
    newBit = 1;
    for column = 1:size(X,2)
        for bit = 1:bitSequence(column)
            bitValue = bitget(X(:,column),bit);
            Y(:,newColumn) = bitset(Y(:,newColumn),newBit,bitValue);
            if newBit == (2^LOG_BYTE_SIZE)
                newColumn = newColumn + 1;
                newBit = 1;
            else
                newBit = newBit + 1;
            end
        end
    end
end