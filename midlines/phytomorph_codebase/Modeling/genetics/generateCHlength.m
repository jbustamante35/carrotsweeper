function [cLv] = generateCHlength(pdfName,para,cN)
    for e = 1:cN
        cLv(e) = round(random(pdfName,para{:},1));
    end
end