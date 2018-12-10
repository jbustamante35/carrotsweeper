D = genD();
%%
n = 3;
for n = 3
    N = n^2;
    tm = clock;
    for e = 1:2^N
        S = myGen3(e);
        S = reshape(fliplr((S(:)')),[3 3]);
        T(e,:) = S(:);
        [P(e) pth(:,:,e)] = vTest3(S,D);
    end
    tmr(n) = etime(clock,tm);
    n
end
%%
[baseFilter f] = generateFilterSet_v2(T,P);
%%
fidx = find(P);
for f = 1:numel(fidx)
    myGen3(fidx(f))
end