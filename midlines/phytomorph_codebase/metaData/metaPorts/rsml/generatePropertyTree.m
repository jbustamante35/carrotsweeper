function [docNode] = generatePropertyTree(docNode,metaNode,tripleCell)

    

    propertySElement = docNode.createElement('property-definitions'); 
    metaNode.appendChild(propertySElement);
    
    
    for e = 1:size(tripleCell,1)
        
        definitionElement = docNode.createElement('property-definition'); 
        
        % label of property
        labelElement = docNode.createElement('label'); 
        value = tripleCell{e,1};
        if ~ischar(value)
            value = num2str(value);
        end
        labelElement.appendChild(docNode.createTextNode(value));
        
        
        % type of property
        typeElement = docNode.createElement('type'); 
        value = tripleCell{e,2};
        if ~ischar(value)
            value = num2str(value);
        end
        typeElement.appendChild(docNode.createTextNode(value));
        
        
        % units of property
        unitsElement = docNode.createElement('unit'); 
        value = tripleCell{e,3};
        if ~ischar(value)
            value = num2str(value);
        end
        unitsElement.appendChild(docNode.createTextNode(value));
        
        
        definitionElement.appendChild(labelElement);
        definitionElement.appendChild(typeElement);
        definitionElement.appendChild(unitsElement);
        
        
        propertySElement.appendChild(definitionElement);
        

    end
end