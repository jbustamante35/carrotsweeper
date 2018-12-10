function [md] = extractMarkerDataFromPop(P)
    md = zeros(size(P(1).ch,1),numel(P));
    for e = 1:numel(P)
        md(:,e) = P(e).ch;
    end
    md = round(md);
end