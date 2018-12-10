function [fileName] = fuse2irods2(fileName,localMountPoint,irodsBase)
    if iscell(fileName)
        for e = 1:numel(fileName)
            fileName{e} = strrep(fileName{e},localMountPoint,irodsBase);        
        end
    elseif ischar(fileName)
        fileName = strrep(fileName,localMountPoint,irodsBase);    
    end
end