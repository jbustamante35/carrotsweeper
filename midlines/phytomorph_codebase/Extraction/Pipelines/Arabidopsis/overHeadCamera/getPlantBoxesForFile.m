function [box] = getPlantBoxesForFile(table,fileName)
    box = getFieldForFileName(fileName,table,'CropBoxes');
    [box] = orderRedCropBoxes(box);
end