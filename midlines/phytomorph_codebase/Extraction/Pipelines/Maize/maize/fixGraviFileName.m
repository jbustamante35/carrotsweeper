function [imagePath] = fixGraviFileName(fileName)
    % load the first image
    [pth,nm,ext] = fileparts(fileName);
    imagePath = strrep(nm,'SPACE',' ');
    imagePath = strrep(imagePath,'SLASH',filesep);
    %imagePath = strrep(imagePath,'-',filesep); 
    %imagePath = strrep(imagePath,[filesep ')'],'-)');
    imagePath = [imagePath filesep];
end