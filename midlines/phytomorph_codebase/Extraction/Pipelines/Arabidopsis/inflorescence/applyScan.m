function [Y] = applyScan(X,net)
    X = permute(X,[1 3 2]);
    sz = size(X);
    X = reshape(X,[prod(sz(1:2)) sz(3)]);
    Y = net.classify(X);
    Y = double(Y)-1;
end