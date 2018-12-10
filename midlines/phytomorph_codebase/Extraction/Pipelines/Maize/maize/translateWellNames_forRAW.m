function [kernelVec genoVec popVec plateVec isoVec indexVec femaleVec maleVec wtsVec] = translateWellNames_forRAW(f)
    import java.util.HashMap;
    kernelVec = HashMap();
    genoVec = HashMap();
    popVec = HashMap();
    plateVec = HashMap();
    isoVec = HashMap();
    indexVec = HashMap();
    femaleVec = HashMap();
    maleVec = HashMap();
    wtsVec = HashMap();
    % constants
    ROWnum = 6;
    COLnum = 8;
    ROWvec = {'A','B','C','D','E','F'};
    COLvec = {'1','2','3','4','5','6','7','8'};
    for e = 1:numel(f.position)
        toB = 0;
        % positions number
        NUMpos = f.position{e};
        % position if out of plate
        if NUMpos > ROWnum*COLnum
            NUMpos = NUMpos - ROWnum*COLnum;
            toB = 1;
        end
        % plate name
        plateN = f.plateName{e};
        plateN = makeKey(plateN);
        if strcmp(plateN,'wisn1013369')
            g = 1;
        end
        fprintf(['Working on plate name:' plateN '\n']);
        try
            if ~isempty(plateN)
                %{
                % A B fix
                if strcmp(plateN(end),'A') | strcmp(plateN(end),'B')  
                    plateN = [plateN(1:end-1) '-' plateN(end)];    
                end
                %}
                
                if toB
                    plateN = [plateN 'b'];
                end
                
                % colN and rowN
                [colN,rowN] = ind2sub([COLnum,ROWnum],NUMpos);
                colV = ['_' COLvec{colN}];
                rowV = ['_' ROWvec{rowN}];

                wellN = [rowV colV];

                
                key1 = [plateN '*' wellN];

                
                % debug - 07s205615
                if strcmp(wellN,'_F_5') & ~isempty(strfind(plateN,'206012'))
                    stop = 1;
                end
                
                % wisn1044335*_E_4
                if strcmp(wellN,'_E_4') & ~isempty(strfind(plateN,'wisn1044335'))
                    stop = 1;
                end
                
                
                
                
                fprintf([key1 '-->' wellN '-->' num2str(NUMpos) '--' num2str(e) '--' num2str(numel(f.position)) '\n']);
                tmpD = kernelVec.get(key1);
                if size(tmpD,2) == 1
                    tmpD = tmpD';
                end
                size(tmpD)
                tmpD = cat(1,tmpD,f.specData(e,:));
                kernelVec.put(key1,tmpD);
                genoVec.put(key1,f.genoType{e});
                popVec.put(key1,f.POP{e});
                plateVec.put(key1,f.plateName{e});
                isoVec.put(key1,f.specIsolate{e});
                indexVec.put(key1,f.specIndex{e});
                femaleVec.put(key1,f.specFemale{e});
                maleVec.put(key1,f.specMale{e});
                
                
                
                tmpD = wtsVec.get(key1);
                if size(tmpD,1) == 1
                    tmpD = tmpD';
                end
                size(tmpD)
                tmpD = cat(1,tmpD,f.specWTS{e});
                wtsVec.put(key1,tmpD);
                
                
            end
        catch ME
           fprintf(['Error on plate name:' plateN '\n']);
        end
    end
end