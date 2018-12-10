function [IDX] = findImageFramePair(metaTable,typeQuery)
    IDX(:,1) = find(strcmp(metaTable.type,typeQuery));
    for e = 1:size(IDX,1)
        keyQuery = metaTable.frameKey(IDX(e,1));
        IDX(e,2) = queryMetaTableKey(metaTable,keyQuery);
    end
end