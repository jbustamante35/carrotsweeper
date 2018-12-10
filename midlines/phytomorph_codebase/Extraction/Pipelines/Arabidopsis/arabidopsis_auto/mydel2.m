function [d] = mydel2(I)
    [g1 g2] = gradient(I);
    [g11 g12] = gradient(g1);
    [g21 g22] = gradient(g2);
    d = g11 + g22;
end