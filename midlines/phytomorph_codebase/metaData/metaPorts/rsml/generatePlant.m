function [docNode] = generatePlant(docNode,plantStruct)
    
    docRootNode = docNode.getDocumentElement;
    sceneElement = docNode.createElement('scene'); 
    docRootNode.appendChild(sceneElement);
    
    for p = 1:numel(plantStruct)
        plantNode = docNode.createElement('plant'); 
        sceneElement.appendChild(plantNode);
        
        curveSet = plantStruct{p}.curveSet;
        functionSets = plantStruct{p}.functionSets;
        
        [docNode] = renderPolyLine(docNode,plantNode,curveSet,functionSets,'root');
    end
end