function [out,pTable] = masterSpore(fileName,patchSZ,oPath)
    out= [];
    [patchStructure,nI] = extractPatches(fileName,patchSZ,100);
    [patchStructure,pTable] = extractPixelListsFromPatches(patchStructure,size(nI),fileName);
    if ~isempty(oPath)
        mkdir(oPath)
        writetable(pTable,[oPath fileName '___phenotypeTable.txt'])
    end
    %writetable(pTable,[oPath fileName '___phenotypeTable.txt'])
    %[out] = paintImage(nI,patchStructure,oPath,fileName);
end