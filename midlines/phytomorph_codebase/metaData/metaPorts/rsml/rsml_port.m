docNode = com.mathworks.xml.XMLUtils.createDocument('rsml');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('xmlns:po','http://www.plantontology.org/xml-dtd/po.dtd');

metaElement = docNode.createElement('metadata'); 
docRootNode.appendChild(metaElement);


versionElement = docNode.createElement('version'); 
versionElement.appendChild(docNode.createTextNode('1'));
metaElement.appendChild(versionElement);


unitElement = docNode.createElement('unit'); 
unitElement.appendChild(docNode.createTextNode('inch'));
metaElement.appendChild(unitElement);


resolutionElement = docNode.createElement('resolution'); 
resolutionElement.appendChild(docNode.createTextNode('30'));
metaElement.appendChild(resolutionElement);
%% si some goodness
for p = 1:3

    numCurves = 3;
    numDims  = 2;
    curveLength = 10;
    for e = 1:numCurves
        tmpCurve = rand(10,numDims);
        curveSet{e} = tmpCurve;

        tmpF = rand(10,1);
        functionSets{e}{1}.name = 'alice';
        functionSets{e}{1}.data = tmpF;

        tmpF = rand(10,1);
        functionSets{e}{2}.name = 'bob';
        functionSets{e}{2}.data = tmpF;
    end
    
    plantS{p}.curveSet = curveSet;
    plantS{p}.functionSets = functionSets;
    
end
%%
docNode = com.mathworks.xml.XMLUtils.createDocument('rsml');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('xmlns:po','http://www.plantontology.org/xml-dtd/po.dtd');
exampleDefsCell{1,1} = 'bob';
exampleDefsCell{1,2} = 'person';
exampleDefsCell{1,3} = 'live';

exampleDefsCell{2,1} = 'alice';
exampleDefsCell{2,2} = 'rock';
exampleDefsCell{2,3} = 'dead';

[docNode] = generateMetaDataTree(docNode,exampleDefsCell);


[docNode] = generatePlant(docNode,plantS);



xmlFileName = ['funtestDone.rsml'];
xmlwrite(xmlFileName,docNode);
type(xmlFileName);

%%



%%
%metaElement.appendChild(metaElement);


xmlFileName = ['funtestDone.xml'];
xmlwrite(xmlFileName,docNode);
type(xmlFileName);




%%

for i=1:20
    thisElement = docNode.createElement('child_node'); 
    thisElement.appendChild... 
        (docNode.createTextNode(sprintf('%i',i)));
    docRootNode.appendChild(thisElement);
end
docNode.appendChild(docNode.createComment('this is a comment'));

xmlFileName = ['funtest.xml'];
xmlwrite(xmlFileName,docNode);
type(xmlFileName);