function [uqX] = generateUQbits(N)
    uqX = [];
    B = uint8(0:(2^N-1))';
    for e = 1:N
        uqX = [uqX;bitget(B,e)'];
    end
    uqX = uqX';
end