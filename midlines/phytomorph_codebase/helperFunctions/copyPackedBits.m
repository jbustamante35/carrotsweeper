function [Y] = copyPackedBits(X,ridx)
    BYTE = floor(ridx/8)+1;
    REM = rem(ridx,8)+1;
    Y = zeros(size(X,1),size(X,2)-floor((1/8)*numel(ridx)),'uint8');
    bitCNT = 1;
    byteCNT = 1;
    for byte = 1:size(X,2)
        for b = 1:8
            if ~any(BYTE==byte .* b==REM)
                Y(:,byteCNT) = bitset(Y(:,byteCNT),bitCNT,bitget(X(:,byte),b));
                if bitCNT == 8
                    bitCNT = 1;
                    byteCNT = byteCNT + 1;
                else
                    bitCNT = bitCNT + 1;
                end
            end
        end
    end
end