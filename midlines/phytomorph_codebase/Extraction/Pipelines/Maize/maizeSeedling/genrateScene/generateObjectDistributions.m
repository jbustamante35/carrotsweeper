function [POS] = generateObjectDistributions(masterTable)
    UQ = unique(masterTable.type);
    %distributions = zeros(sceneSize,numel(UQ),2);
    
    for u =1:numel(UQ)
        X = [];
        SIZE = [];
        fidx = find(strcmp(masterTable.type,UQ{u}));
        
        X(:,2) = masterTable.boundingBox1(fidx);
        X(:,1) = masterTable.boundingBox2(fidx);
        
        
        SIZE(:,2) = masterTable.boundingBox3(fidx);
        SIZE(:,1) = masterTable.boundingBox4(fidx);
        
        POS.CM{u} = X + SIZE/2;
        POS.LOC{u} = X;
        POS.LABEL{u} = UQ{u};
        
        
    end
end