function [] = obtainFeatureMaps(rawData,featureBank,keyList,N,oPort)
    parfor e = 1:N
        
        % register toCompute
        toCompute = [];
        [fileKey] = getFileKey(rawData,e);
        for key = 1:numel(keyList)
            outFile = [oPort fileKey '--' keyList{key} '.mat'];
            toCompute(key) = ~exist(outFile);
        end
        
        if any(toCompute)
            % get the data
            [data fileKey] = getData(rawData,e);
            % compute if needed
            for key = 1:numel(keyList)
                outFile = [oPort fileKey '--' keyList{key} '.mat'];
                if ~exist(outFile);
                    func = featureBank.getFunction(keyList{key});
                    fm = func(data,e);
                    
                    storeFeatures(fileKey,keyList{key},fm,oPort);
                end
            end
        end
    end
end

function [data fileKey] = getData(rawData,e)
    tim = clock;
    inFile = rawData{e}.fileName;
    fprintf(['Starting load of ' inFile  '\n']);
    [pth fileKey ext] = fileparts(inFile);
    data = myReader(rawData{e});
    fprintf(['Ending load of ' inFile '\n']);
end

function [fileKey] = getFileKey(rawData,e)
    inFile = rawData{e}.fileName;
    [pth fileKey ext] = fileparts(inFile);
end

function [] = storeFeatures(fileKey,featureKey,featureMap,oPort)
    tim = clock;
    outFile = [oPort fileKey '--' featureKey '.mat'];
    fprintf(['Starting save of ' outFile '\n']);
    data = featureMap.fV;
    index = featureMap.index;
    sz = featureMap.sz;
    save(outFile,'data','index','sz');
    fprintf(['Ending save of ' outFile '\n']);
end
