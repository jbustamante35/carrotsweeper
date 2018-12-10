function [c] = particleCountingGauge(molecule,c)
    if isempty(c)
        c = 0;
    end
    c = c + size(molecule,2);
end