function [out] = jarrayTocell(in)
    out = {};
    for e = 1:in.size()
        out{e} = in.get(e-1);
    end
end