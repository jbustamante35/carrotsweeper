function [snip] = extractSnip(S,NP)
    snip = extractSnipOverStack(S,NP);
end

function [snip] = extractSnipOverStack(S,NP)
    for e = 1:numel(S)
        snip(:,:,e) = extractSnipOverObjects(S(e).midlines.data,NP);
        
    end
end
    