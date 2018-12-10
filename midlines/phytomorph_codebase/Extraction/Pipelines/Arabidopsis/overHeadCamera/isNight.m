function [tf] = isNight(I,thresh)
    tf = mean(I(:)) < thresh;
end