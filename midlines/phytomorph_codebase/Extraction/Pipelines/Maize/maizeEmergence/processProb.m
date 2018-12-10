function [out] = processProb(in,label)
    oin = in;
    %fin = in > .95;
    fin = in > .9*max(in);
    fin = imfilter(double(fin),fspecial('gaussian',[1 20],7));
    in = in.*fin;
    sin = imfilter(in,fspecial('average',[1 20]));
    ds = gradient(sin);
    [~,pidx] = max(ds);
    out = zeros(size(in));
    out(pidx:end) = 1;
    
    
    fidx = find(label);
    pidx = find(out);
    if ~isempty(pidx) & ~isempty(fidx)
        delta = abs(pidx(1) - fidx(1));
        if delta > 10
            stop = 1;
        end
    end
end