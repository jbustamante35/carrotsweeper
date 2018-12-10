function [T] = meltTensor(T)
    sz = size(T);
    T = reshape(T,[prod(sz(1:end-1)) sz(end)])';
end