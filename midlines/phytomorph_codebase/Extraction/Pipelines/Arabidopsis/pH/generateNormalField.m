function [NOR] = generateNormalField(dB,segmentSize)
    NOR = [];
    segmentSize = 31;
    tmp = dB{1}';
    for e1 = 1:size(tmp,2)
        % get curve segment
        sig = tmp(:,1:segmentSize)';
        % fit with least squares spline
        fn = spap2(1,2,[1:size(sig,1)]',sig');
        % get function derivative
        fn2 = fnder(fn,1);
        pvec = fnval(fn,(segmentSize-1)/2);
        % eval the der for tangent vec
        tvec = fnval(fn2,(segmentSize-1)/2);
        tvec = tvec/norm(tvec);
        nvec = [tvec(2);-tvec(1)];
        tmp = circshift(tmp,[0 -1]);
        NOR = [NOR -nvec];
    end
    NOR = circshift(NOR,[0 (segmentSize-1)/2])';
end