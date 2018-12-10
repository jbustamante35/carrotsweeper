function [out] = parseInput(inputString)
    delimiter = '**';
    fidx = strfind(inputString,'**');
    offSet = numel(delimiter);
    type = inputString(fidx(1)+offSet:fidx(2)-1);
    value = inputString(fidx(2)+offSet:fidx(3)-1);
    if ~isempty(strfind(value,'matfile'))
        fidx = strfind(value,',');
        matfileName = value(1:fidx(1)-1);
        matfileName(1:8) = [];
        varsavedName = value(fidx(1)+1:end);
        data = load(matfileName,varsavedName);
        out = data.(varsavedName);            
    end
end
