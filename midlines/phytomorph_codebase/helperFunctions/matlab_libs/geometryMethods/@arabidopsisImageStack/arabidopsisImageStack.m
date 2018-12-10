classdef arabidopsisImageStack 
    properties
        imageStack;
        count;
    end
    
    methods
        
        function [obj] = arabidopsisImageStack(imageStack)
            for e =1:numel(imageStack)
                obj.imageStack{e} = arabidopsisImage(imageStack{e});
            end
        end
        
        function [] = flush(obj)
            for e = 1:numel(obj.imageStack)
                obj.imageStack{e}.imageData = [];
            end
        end
        
        function [name] = getName(obj)
            [p n e] = fileparts(obj.imageStack{1}.filename);
            name = strrep(p,'/','--');
        end
            
        function [curveSequence] = extractContourSequences(obj)
            % extract the kernel centers
            kernelCenters = obj.imageStack{1}.extractKernelCenterPoints();
            % extract the initial contours
            curveSet = obj.imageStack{1}.extractClosedContours(@(x)maizeContour(x),10,1,[100 500]);
            % extract full curve set from the 
            curveSet = curveStore.containsPoints(curveSet,kernelCenters,@(x)any(x));
            % generate containment map for curveSet
            containmentMap = curveStore.containmentMap(curveSet);
            % for each "root (graph def)" object generate a representation 
            repCurveSet = curveStore.generateRepCurve(containmentMap,curveSet,5);
            % sort the kernels by the kernel centers
            for k = 1:size(kernelCenters,2)
                for e = 1:numel(repCurveSet)
                    if ~isempty(curveStore.containsPoints(repCurveSet(e),kernelCenters(:,k),@(x)any(x)))
                        sidx(k) = e;
                    end
                end
            end
            repCurveSet = repCurveSet(sidx);
            
            
            for e = 1:numel(repCurveSet)
                curveSequence{e} = maizeContourSequence(obj);
                curveSequence{e}.putContour(repCurveSet(e));
                for f = 1:numel(obj.imageStack)-1
                    tmpSequence = curveSequence{e}.getContour(f);
                    BOX = tmpSequence.getPadBoundingBox(50);
                    OFFSET = BOX(1:2)';
                    I = obj.imageStack{f+1}.readCropped(BOX);                    
                    tmpContours = I.extractClosedContours(@(x)maizeContour(x),10,1,[100 500]);
                    tmpContours = curveStore.containsPoints(tmpContours,kernelCenters(:,e) - OFFSET,@(x)any(x));                                       
                    % match the curve                    
                    dist = [];
                    for c = 1:numel(tmpContours)
                        dist(c) = tmpSequence.meanDistance(tmpContours(c));
                    end
                    [J midx] = min(dist);
                    tmpContours(midx(1)).displace(OFFSET);
                    curveSequence{e}.putContour(tmpContours(midx(1)));
                    %{
                    
                    h = imshow(I.imageData,[]);                    
                    hold on
                    curveSequence{e}(f+1).plot(h);
                    drawnow
                    %}
                    %curveSequence{e}(f+1).displace(OFFSET);
                    fprintf(['done@curve:' num2str(e) ':image:' num2str(f) '\n'])
                end
            end
        end
    end
end