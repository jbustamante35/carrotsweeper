function [in] = removeDoubleQuotes(in)
    if isdeployed()        
        in(1) = [];
        in(end) = [];
    end
end