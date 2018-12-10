function [structList] = appendStruct(structList,newStruct)
    if isempty(structList)
        clear structList;
        structList(1) = newStruct;
    else
        Nflds = fieldnames(newStruct);
        defaultValue = NaN;
        for s = 1:numel(structList)
            flds = fieldnames(structList(s));

            toAddtoNew = setdiff(flds,Nflds);
            toAddtoOld = setdiff(Nflds,flds);

            cDoc = structList(s);

            for f = 1:numel(toAddtoOld)
                tmp = newStruct.(toAddtoOld{f});
                sz = tmp.size;
                type = tmp.type;
                if strcmp(type,'double')
                    tmpD = NaN*zeros(sz);
                else
                    tmpD = 'NaN';
                end

                cDoc = generatePhenotypeNode(cDoc,tmpD,tmp.name,toAddtoOld{f});
            end


            for f = 1:numel(toAddtoNew)
                tmp = cDoc.(toAddtoNew{f});

                sz = tmp.size;
                type = tmp.type;
                if strcmp(type,'double')
                    tmpD = NaN*zeros(sz);
                else
                    tmpD = 'NaN';
                end
                newStruct = generatePhenotypeNode(newStruct,tmpD,tmp.name,toAddtoNew{f});
            end


            newList(s) = cDoc;















        end
        newList(end+1) = newStruct;
        structList = newList;
    end
end