function [partP] = genFullPath(partP,mainP)    
    % attach main path to the file name
    mainP = attachSlash(mainP);
    for i = 1:size(partP,1)
        partP(i).name = [mainP partP(i).name];
    end
end