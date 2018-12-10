function [GMModel] = makeGMMset(S)
    sz = size(S);
    S = permute(S,[3 1 2 4]);
    sz = size(S);
    S = reshape(S,[prod(sz(1)) prod(sz(2:4))]);
    options = statset('Display','iter');
    GMModel = fitgmdist(S',4,'Start','plus','Options',options);
end