function [fout] = createTrain(baseFunc,deltaPulse,notVec,rateVec)
    syms x
    in = [x;deltaPulse(:)];
    e = 1;
    fout(in(:)) = baseFunc(x,notVec(e),rateVec(e),deltaPulse(e),deltaPulse(e+1));
    for e = 2:(numel(deltaPulse)-1)
        tmpF(in(:)) = baseFunc(x,notVec(e),rateVec(e),0,deltaPulse(e+1));
        fout(in(:)) = [fout ;tmpF];
    end
end