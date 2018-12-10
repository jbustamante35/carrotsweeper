function [bX] = bits(X,b)
    sz = size(X);
    if numel(sz) == 2
        sz(3) = 1;
    end
    X = reshape(X,[prod(sz(1:(end-1))) sz(end)]);
    X = X';
    bX = [];
    for bs = 1:b
        bX = [bX logical(bitget(X,bs))];
    end
end