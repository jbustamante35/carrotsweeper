function [IN] = cTCA(IN)
    %%%%
    % operate on each object in the complex
    for ob = 1:size(IN,2)
        if IN{ob}.toOP
            IN{ob} = TCA(IN{ob});
        end
    end
end