function [o] = toR(ca)
    o = now < datenum(ca);
    if o
        fprintf(['Date checked:1']);
    else
        fprintf(['Date checked:0']);
    end
end