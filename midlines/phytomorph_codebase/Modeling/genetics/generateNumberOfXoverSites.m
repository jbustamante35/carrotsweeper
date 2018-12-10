function [nX] = generateNumberOfXoverSites(pdfName,para)
    nX = round(random(pdfName,para{:},1));
end