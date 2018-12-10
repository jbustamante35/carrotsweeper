function [T,NM] = stomataExtractionMolecules(patchData,extractType,extractArgs)
    switch extractType
        case 'stomata_fft'
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extractArgs{1} = NF - number of frequencies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M = mfftm(size(patchData,2),extractArgs{1});
            patchData = permute(patchData,[2 1]);
            fF = mtimesx(M,patchData);
            fF = permute(fF,[2 1]);
            T{1} = abs(fF);
            T{2} = angle(fF);
            
            
            T = freezeTensor(T,false);
            NM = 'T';
    case 'stomata_fft_compress'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extractArgs{1} = NF - number of frequencies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M = mfftm(size(patchData,2),extractArgs{1});
            patchData = permute(patchData,[2 1]);
            fF = mtimesx(M,patchData);
            fF = permute(fF,[2 1]);
            T{1} = abs(fF);
            T{2} = angle(fF);
            data.T = freezeTensor(T,false);
            data = prepareData(data,extractArgs{2},extractArgs{3});
            T = data.C;
            NM = 'C';
        case 'fullApply'
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extractArgs{1} = NF - number of frequencies
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            M = mfftm(size(patchData,2),extractArgs{1});
            patchData = permute(patchData,[2 1]);
            fF = mtimesx(M,patchData);
            fF = permute(fF,[2 1]);
            T{1} = abs(fF);
            T{2} = angle(fF);
            data.T = freezeTensor(T,false);
            data = prepareData(data,extractArgs{2},extractArgs{3});
            [probMap] = applyAIlayer(data,extractArgs{4},extractArgs{2},extractArgs{3},[1 1]);
            T = squeeze(probMap);
            NM = 'P';
    end
end