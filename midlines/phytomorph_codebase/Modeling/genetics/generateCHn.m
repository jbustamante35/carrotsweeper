function [cN] = generateCHn(pdfName,para,N)
    cN = round(random(pdfName,para{:},N{:}));
end