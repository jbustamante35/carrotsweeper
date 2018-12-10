function [imageStack] = removeANDsort(imageStack)
    %% remove not compliant images - those that are not numbers
    n = [];
    for e = 1:numel(imageStack)
        try
            [p,tn,ex] = fileparts(imageStack{e});
            rm(e) = 0;
            if strcmp(tn(1),'.')
                rm(e) = 1;
            end
        catch
            rm(e) = 1;
        end
    end
    imageStack(find(rm)) = [];
    %% sort the images
    n = [];
    for e = 1:numel(imageStack)
        try
            [p,tn,ex] = fileparts(imageStack{e});
            n(e) = str2num(tn);    

        catch

        end
    end
    [~,sidx] = sort(n);
    imageStack = imageStack(sidx);
end