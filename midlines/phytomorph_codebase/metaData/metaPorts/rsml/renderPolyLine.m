function [docNode] = renderPolyLine(docNode,plantNode,curveSet,functionSets,polyLineType)
    % polyLineType := root is default - but this is tops
    
    
    
    
    
    for e = 1:numel(curveSet)
        % create organ
        organElement = docNode.createElement(polyLineType);
        % set the ID as unique
        organElement.setAttribute('ID',num2str(e));
        % attach organ to plant
        plantNode.appendChild(organElement);
        % attach geometry to organ
        geometryElement = docNode.createElement('geometry'); 
        organElement.appendChild(geometryElement);
        
       
        
        % attach polyline to geometry
        polylineElement = docNode.createElement('polyline'); 
        geometryElement.appendChild(polylineElement);
        % create the points inthe poly line
        for pt = 1:size(curveSet{e},1)
            pointElement = docNode.createElement('point');
            xValue = num2str(curveSet{e}(pt,1));
            pointElement.setAttribute('x',xValue);
            yValue = num2str(curveSet{e}(pt,2));
            pointElement.setAttribute('y',yValue);
            if size(curveSet{e},2) == 3
                zValue = num2str(curveSet{e}(pt,3));
                pointElement.setAttribute('z',zValue);
            end
            polylineElement.appendChild(pointElement);
        end
        
        % attach any functions
        functionsNode = docNode.createElement('functions'); 
        organElement.appendChild(functionsNode);
        % on fths function of total
        for f = 1:numel(functionSets{e})
            tmpFunctionNode = docNode.createElement('function'); 
            tmpFunctionNode.setAttribute('name',functionSets{e}{f}.name);
            tmpFunctionNode.setAttribute('domain','polyline');
            % pth point along the eth curve and the fth function
            for pt = 1:size(functionSets{e}{f}.data,1)
                 tmpSample = docNode.createElement('sample');
                 value = functionSets{e}{f}.data(pt,:);
                 if ~ischar(value)
                    value = num2str(value);
                    value = strrep(value,'     ',',');
                 end
                 tmpSample.appendChild(docNode.createTextNode(value));
                 tmpFunctionNode.appendChild(tmpSample);
            end
            functionsNode.appendChild(tmpFunctionNode);
        end
        
    end
end