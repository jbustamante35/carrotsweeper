function [res] = accumulatorFunc(particleSpray,func,res)
    if isempty(res)
        res = 0;
    end
    res = res + func(particleSpray);
end