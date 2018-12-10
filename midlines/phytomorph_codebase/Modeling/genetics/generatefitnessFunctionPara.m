function [fitness] = generatefitnessFunctionPara(pdfName,cLv,maxSites,maxMarker)
    % pdfName := distribution to pull from
    % cLv := chromosome length vector
    % maxSites := max number of sites to select from
    for e = 1:numel(cLv)
        para = {1 cLv(e)};
        fitness(e).numSites = randi(maxSites,1,1);
        fitness(e).sites = round(random(pdfName,para{:},fitness(e).numSites,1));
        fitness(e).values = randi(maxMarker,fitness(e).numSites,1);
    end
end