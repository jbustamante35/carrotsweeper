classdef phG < handle
    properties
        phNl;
        phEl;
        pointTable = [];
        edgeTable;
    end
    
    methods
        function [obj] = phG(points)
            obj.pointTable = points;
            obj.phNl = phN.empty();
            obj.phEl = phE.empty();
            obj.edgeTable = sparse(size(points,1),size(points,1));
        end
        
        function [obj] = addNode(obj,n)
            obj.phNl(numel(obj.phNl)+1) = n;
        end
        
        
        function [obj] = addEdge(obj,e)
            obj.phEl(numel(obj.phEl)+1) = e;
            obj.edgeTable(e.nodeL(1),e.nodeL(2)) = obj.edgeTable(e.nodeL(1),e.nodeL(2)) + 1;
            obj.edgeTable(e.nodeL(2),e.nodeL(1)) = obj.edgeTable(e.nodeL(2),e.nodeL(1)) + 1;
        end
        
        
        function []  = displayNodes(obj,CL,idx)
            hold on
            if nargin < 3
                idx = 1:numel(obj.phNl);
            end
            if nargin <= 1
                CL = 'go';
            end
            for e = 1:numel(idx)
                plot(obj.pointTable(obj.phNl(idx(e)).index,1),obj.pointTable(obj.phNl(idx(e)).index,2),CL)
            end
        end
        
        function [] = displayEdges(obj)
            for e = 1:numel(obj.phEl)
                vec = obj.pointTable(obj.phEl(e).nodeL,:);
                plot(vec(:,1),vec(:,2),'r');
            end
        end
        
        function [fidx] = findEndPoints(obj)
            fidx = find(sum(obj.edgeTable,1) == 1);
        end
        
        % find non end points
        function [fidx] = findnEndPoints(obj)
             fidx = find(sum(obj.edgeTable,1) >= 2);
        end
        
        function [] = displayEndPoints(obj)
            fidx = obj.findEndPoints;
            hold on
            for e = 1:numel(fidx)
                plot(obj.pointTable(fidx(e),1),obj.pointTable(fidx(e),2),'bo')
            end
        end
        
        function [t] = traceFrom(obj,s,edgeTable)
            if nargin <=2
                edgeTable = obj.edgeTable;
            end
            t = edgeTable*s;
        end
        
        function [startPoint] = traceToEnd(obj,startPoint,maxS)
            if nargin <= 2
                maxS = inf;
            end
            edgeTable = obj.edgeTable;
            fidx = obj.findEndPoints();
            % create self links for end points
            for e = 1:numel(fidx)
                edgeTable(fidx(e),fidx(e)) = 1;
            end
            flag = 1;
            e = 1;
            disp = 0;
            while flag
                startPoint(:,e+1) = obj.traceFrom(startPoint(:,e),edgeTable);
                % find forward flow to prevent backflow
                ep = find(startPoint(:,e+1));
                sp = find(startPoint(:,e));
                for e1 = 1:numel(ep)
                    for e2 = 1:numel(sp)
                        if sp(e2) ~= ep(e1)
                            edgeTable(sp(e2),ep(e1)) = 0;
                        end
                    end
                end
                %startPoint(any(sum(startPoint,2) >= 2,2),e+1) = 0;
                
                if disp
                    close all
                    hold on
                    obj.displayEdges();
                    obj.displayEndPoints();
                    obj.displayNodes();
                    cur = find(startPoint(:,e));
                    plot(obj.pointTable(cur,1),obj.pointTable(cur,2),'k*');
                    cur = find(startPoint(:,e+1));
                    plot(obj.pointTable(cur,1),obj.pointTable(cur,2),'c*');
                    hold off
                    drawnow
                end
                if all(startPoint(:,e+1) == startPoint(:,e)) | maxS == e;
                    flag = 0;
                end
             
                e = e + 1;
                %waitforbuttonpress
            end
            rm = all(startPoint==0,1);
            startPoint(:,rm) = [];
        end
    end
end