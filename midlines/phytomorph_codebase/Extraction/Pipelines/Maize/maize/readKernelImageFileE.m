function [I] = readKernelImageFileE(fileName,ext)
    fileName = strrep(fileName,'/_','/');
    I = tryOne(fileName,ext);
    if ~I
        fileName = strrep(fileName,'/mnt/spaldingdata/','/mnt/snapper/');
        I = tryOne(fileName,ext);
        if ~I
            fprintf(['Not Found:' fileName '\n']);
        else
            fprintf(['Found:' fileName '\n']);
        end
    else
        fprintf(['Found:' fileName '\n']);
    end
    %{
    if ~isempty(I)
        fileName
    else

    end
    %}
end

function [tmp] = tryOne(fileName,ext)
    tmp = [];
    try
        tmp = exist([fileName ext]);
    if ~tmp
        
            rawName = fileName;
            if strcmp(ext,'_S2.tiff')
                rawName = strrep(rawName,'_F.tiff',ext);
            end
            tmp = exist(rawName);
        if ~tmp
            
                tmp = exist([strrep(fileName,'gray/07S','gray/07') ext]);
            if ~tmp
                 
                    rawName = [strrep(fileName,'gray/07S','gray/07')];
                    if strcmp(ext,'_F.tiff')
                        rawName = strrep(rawName,'_S2.tiff',ext);
                    end
                    tmp = exist(rawName);
                 if ~tmp
                    
                        rawName = [strrep(fileName,'_WIS','WIS')];
                        if strcmp(ext,'_F.tiff')
                            rawName = strrep(rawName,'_S2.tiff',ext);
                        end
                        tmp = exist([rawName ext]);
                    
                 end
            end
        end
    end
    catch
    end
end