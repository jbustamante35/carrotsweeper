function [oxs xVec] = generateXoverLocs(pdfName,N,para)
    try
        oxs = round(random(pdfName,para{:},N{:}));
        oxs = sort(oxs);
        xVec = zeros(para{end},1)';
        if ~isempty(oxs)
            xVec(oxs) = 1;
        end
        xVec = [xVec xVec];
    catch ME
        ME
    end
end