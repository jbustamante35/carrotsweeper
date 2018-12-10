function [idx] = queryMetaTableKey(metaTable,keyIndex)
    idx = find(metaTable.key == keyIndex);
end