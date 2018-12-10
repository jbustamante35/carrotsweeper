function [M] = collateJSON(FileList,oPath)



    % get list of fields
    data = loadjson(FileList{1});
    docField = fields(data(1));
    baseFields = fields(data.(docField{1}){1});
    lookUps = {};
    sz = [];
    for f = 1:numel(baseFields)
        if ~ischar(data.(docField{1}){1}.(baseFields{f}))
            tmp = data.(docField{1}){1}.(baseFields{f});
            sz = [sz numel(tmp)];
            lookUps{end+1} = baseFields{f};
        end
    end
    
    
    
    % populate the headers
    M = {};
    cnt = 1;
    M{1,1} = 'QR_data';
    M{1,2} = 'Object_number';
    colSTR = 2;
    for f = 1:numel(lookUps)
        for s = 1:sz(f)
            M{1,colSTR+cnt} = lookUps{f};
            cnt = cnt + 1;
        end
    end
    
    dcnt = 1;
    for e = 1:numel(FileList)
        fprintf(['starting document:' num2str(e) ':' num2str(numel(FileList)) '\n']);tic
        [pth,nm,ext] = fileparts(FileList{e});
        fidx = strfind(nm,'_');
        fName = nm(1:(fidx(1)-1));
        
        data = loadjson(FileList{e});
        numDocs = numel(data.(docField{1}));
        for o = 1:numDocs
            cnt = 1;
            M{1+dcnt,1} = fName;
            M{1+dcnt,2} = num2str(o);
            for f = 1:numel(lookUps)
                d = data.(docField{1}){o}.(lookUps{f});
                for i = 1:numel(d)
                    M{1+dcnt,colSTR+cnt} = d(i);
                    cnt = cnt + 1;
                end
            end
            dcnt = dcnt + 1;
        end
        fprintf(['ending document:' num2str(e) ':' num2str(numel(FileList)) ':' num2str(toc) '\n'])
    end
    
    
    
    mkdir(oPath)
    fileOut = [oPath 'masterSheet.csv'];
    cell2csv(fileOut,M);
end