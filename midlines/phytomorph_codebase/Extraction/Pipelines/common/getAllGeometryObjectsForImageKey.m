function [IDX] = getAllGeometryObjectsForImageKey(metaTable,keyIndex)
    % get the table index
    idx = find(metaTable.key == keyIndex);
    % get the frame for the keyIndex
    keyQuery = metaTable.frameKey(idx);
    % find all objects that are in this images frame
    allIDX = find(metaTable.frameKey == keyQuery);







    rm = contains(lower(metaTable.type(allIDX)),'image');
    allIDX(rm) = [];
    metaTable(allIDX,:)


end