function [uKey] = genUQkey(fileName,sChar,nDepth)
    % parse input for "flags" for unique names
    %[p,n,e] = fileparts(fileName);
    
    p = fileName;
    if ~strcmp(p(end),filesep)
        p = [p filesep];
    end
    % find the slash operator
    fidx = strfind(p,filesep);
    
    fidx(1:4) = [];
    
    % generate new base
    new_basePath = p(fidx(1)+1:end);
    % key 
    uKey = strrep(new_basePath,filesep,sChar);
end