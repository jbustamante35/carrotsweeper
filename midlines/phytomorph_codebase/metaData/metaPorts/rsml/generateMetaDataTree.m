function [docNode] = generateMetaDataTree(docNode,tripleCell)
    defFile = 'https://de.cyverse.org/dl/d/9B91681E-AC78-43E9-8ED6-ACC77B49664C/rsml_metadata.csv';
    localFile = [tempname '.csv'];
    websave(localFile,defFile);
    defs = readtext(localFile);
    
    docRootNode = docNode.getDocumentElement;
    
    
    metaElement = docNode.createElement('metadata'); 
    docRootNode.appendChild(metaElement);
    
    for e = 1:size(defs,1)
        tmpElement = docNode.createElement(defs{e,1});
        value = defs{e,2};
        if ~ischar(value)
            value = num2str(value);
        end
        tmpElement.appendChild(docNode.createTextNode(value));
        metaElement.appendChild(tmpElement);
    end
    
    [docNode] = generatePropertyTree(docNode,metaElement,tripleCell);


end