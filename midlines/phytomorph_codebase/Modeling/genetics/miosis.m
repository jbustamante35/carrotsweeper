function [G] = miosis(xM,C)
    G.ch = xM*C.ch;
    G.cn = C.cn;
end