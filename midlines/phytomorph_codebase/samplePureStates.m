function [INT] = samplePureStates(MAX_B,TALL)
    INT = double(randi((2^MAX_B-1),TALL,1)-1);  
end