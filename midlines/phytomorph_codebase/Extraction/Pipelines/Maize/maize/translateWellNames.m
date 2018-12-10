function [kernelVec genoVec popVec kernelID rkernelID] = translateWellNames(f)
    import java.util.HashMap;
    kernelVec = HashMap();
    genoVec = HashMap();
    popVec = HashMap();
    kernelID = HashMap();
    rkernelID = HashMap();
    for e = 1:numel(f.position)
       
        ROWnum = 6;
        COLnum = 8;
        ROWvec = {'A','B','C','D','E','F'};
        COLvec = {'1','2','3','4','5','6','7','8'};
        
    
    
        NUMpos = f.position{e};
        
        
        if NUMpos > ROWnum*COLnum
            NUMpos = NUMpos - ROWnum*COLnum;
        end
        plateN = f.plateName{e};
       
        if strcmp(plateN(end),'A') | strcmp(plateN(end),'B')  
            plateN = [plateN(1:end-1) '-' plateN(end)];    
        end
        
        
        %[rowN,colN] = ind2sub([ROWnum,COLnum],NUMpos);
        [colN,rowN] = ind2sub([COLnum,ROWnum],NUMpos);
        %colN = rem(NUMpos,COLnum)+1;
        colN
        colV = ['_' COLvec{colN}];

        %rowN = floor(NUMpos/COLnum)+1
        rowN
        rowV = ['_' ROWvec{rowN}];
        
        
    
    
    
        wellN = [rowV colV];

        key1 = [plateN '*' wellN];
    
        fprintf([key1 '-->' wellN '-->' num2str(NUMpos) '\n']);
        %kernelVec.put(key1,[f.specData(e,:) f.prediction(e,:)]);
        kernelVec.put(key1,[f.specData(e,:)]);
        
        kernelID.put([key1 '--' f.POP{e}],f.kernel_id{e});
        
        
        rkernelID.put(f.kernel_id{e},[key1 '--' f.POP{e}]);
        genoVec.put(key1,f.genoType{e});
        popVec.put(key1,f.POP{e});
    end
end