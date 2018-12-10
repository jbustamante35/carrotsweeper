function [result] = applyNetworkToPoint(I,P,preprocessFunc,networkToApply,dataBasis,dataMean)
    data = preprocessFunc(I,P);
    data = permute(data,[2 1 3]);
    if nargin > 4
        dataNew = [];
        for d = 1:size(data,3)
            dataC = PCA_REPROJ_T(data(:,:,d),dataBasis{d},dataMean{d});
            sim = PCA_BKPROJ_T(dataC,dataBasis{d},dataMean{d});
            err = sum((data(:,:,d) - sim).^2,1).^.5;
            dataNew = [dataNew;dataC;err];
        end
        data = dataNew;
    end
    result = networkToApply.predict({data});
end