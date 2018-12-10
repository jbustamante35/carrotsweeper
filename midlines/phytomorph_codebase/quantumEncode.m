function [X] = quantumEncode(MAX_B,INT)
    % 0 - 255
    
    BYTE_SIZE = 3;                              % log size of a byte
    MAX_WID = max(1,2^MAX_B / 2^BYTE_SIZE);     % number of bytes to allocate wide


    TALL = size(INT,1);


    X = zeros(TALL,MAX_WID,'uint8');

    BYTE = floor(double(INT)*(2^BYTE_SIZE)^-1)+1;
    %BYTE = rem(double(INT),(2^MAX_B)+1)+1;
    REM = rem(INT,2^BYTE_SIZE)+1;
    POS = sub2ind(size(X),(1:size(X,1))',BYTE);
    X(POS) = bitset(X(POS),REM);
end