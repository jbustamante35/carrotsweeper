function [I] = readKernelImageFile(fileName,ext)
    fileName = strrep(fileName,'/_','/');
    I = tryOne(fileName,ext);
    if isempty(I)
        fileName = strrep(fileName,'/mnt/spaldingdata/','/mnt/snapper/');
        I = tryOne(fileName,ext);
        if ~isempty(I)
            snapper = 1
        end
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
        tmp = imread([fileName ext]);
    catch ME2
        try
            rawName = fileName;
            if strcmp(ext,'_S2.tiff')
                rawName = strrep(rawName,'_F.tiff',ext);
            end
            tmp = imread(rawName);
        catch ME2
            try
                tmp = imread([strrep(fileName,'gray/07S','gray/07') ext]);
            catch ME2
                 try
                    rawName = [strrep(fileName,'gray/07S','gray/07')];
                    if strcmp(ext,'_F.tiff')
                        rawName = strrep(rawName,'_S2.tiff',ext);
                    end
                    tmp = imread(rawName);
                 catch ME2
                    try
                        rawName = [strrep(fileName,'_WIS','WIS')];
                        if strcmp(ext,'_F.tiff')
                            rawName = strrep(rawName,'_S2.tiff',ext);
                        end
                        tmp = imread([rawName ext]);
                    catch ME3

                    end
                 end
            end
        end
    end
end