function [ dim ] = fractalDim( tBin )
%FRACTALDIM Calculates FD using box counting method

    [d1 d2] = size(tBin);
    maxDim = max(size(tBin));
    pad1 = 2^nextpow2(maxDim) - d1;
    pad2 = 2^nextpow2(maxDim) - d2;

    tPad = padarray(tBin, [pad1 pad2], 'post');
    [d1 d2] = size(tPad);
    
    nBox = (size(tPad, 1) * size(tPad, 2)) / (d1*d2);
    
    n = [];
    s = [];
    i = 1;
    while mod(d1, 1) == 0 && mod(d2, 1) == 0
        nBox = (size(tPad, 1) * size(tPad, 2)) / (d1*d2);
        n(i)=0;
        s(i)= sqrt(nBox);
        for r = 0:(sqrt(nBox)-1)
            for c = 0:(sqrt(nBox)-1)
                box = tPad((d1*r + 1):(d1*(r+1)), (d2*c + 1):(d2*(c+1)));
                count = sum(sum(box));
           
                if count > 0
                    n(i) = n(i) +1; 
                end
            end
        end
        
        d1 = d1 / 2;
        d2 = d2 / 2;
        i = i + 1;
    end

    reg = fit(log(s'), log(n'), 'poly1');
    dim = reg.p1;
    
end

