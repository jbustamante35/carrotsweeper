classdef maizeContourSequence < handle
    
    properties
        contourSequence;
        imageStack;
    end
    
    
    methods
        function [obj] = maizeContourSequence(imageStack)
            obj.imageStack = imageStack;
        end
        
        function [ret] = getContour(obj,varargin)
            if nargin == 1;index = numel(obj.contourSequence);
            else;index = varargin{1};end
            ret = obj.contourSequence{index};
        end
        
        function [] = putContour(obj,varargin)
            if nargin == 2;index = numel(obj.contourSequence)+1;
            else;index = varargin{2};end
            obj.contourSequence{index} = varargin{1};
        end
        
        function [] = persist(obj,sName)
            save([sName],'obj');
        end
        
        function [name] = getName(obj)
            [p n e] = fileparts(obj.imageStack.imageStack{1}.filename);
            name = strrep(p,'/','--');
        end
        
        function [] = plot(obj)
            for e = 1:numel(obj.contourSequence)                
                [I offset] = obj.contourSequence{e}.getImage();
                obj.contourSequence{e}.displace(-offset);
                imshow(I,[]);
                hold on;
                obj.contourSequence{e}.plot();
                hold off;
                drawnow;
                obj.contourSequence{e}.displace(offset);
                %waitforbuttonpress
            end
        end
        
        function [curveAtoms] = atomizeSingleContour(obj,segmentSize,index)
            [curveAtoms] = obj.contourSequence{index}.atomizeCurve(segmentSize,@()closedCurveSegment());
        end
        
        function [] = traceMidline(obj,segmentSize,E,U)
            try
                for tm = 1:numel(obj.contourSequence)
                    % atomize the contour
                    atomicFragments = obj.atomizeSingleContour(segmentSize,tm);
                    % stack the sub-atomic fragments
                    tmp = [];
                    for e = 1:numel(atomicFragments)
                        tmp = [tmp;atomicFragments(e).Osegment(:)'];
                    end
                    % project into eigenspace
                    tC = PCA_REPROJ(tmp,E,U);
                    % clear the fragments that are less than the center kernel
                    tC = tC(:,1).*(obj.contourSequence{tm}.centerVec(1,:)  < 0)';
                    [J fidx] = max(tC);
                    % tag the tip
                    obj.contourSequence{tm}.tagTip(fidx);


                    if tm == 1
                        ip = obj.contourSequence{1}.getIntersectionProfile();
                        ip = sum(ip.*ip,1).^.5;
                        d1 = cwt(ip,[30],'gaus1');
                        d1 = abs(d1);
                        jidx = nonmaxsuppts(d1, 10);
                        % threshold on normal-distance metric
                        jidx = jidx .* (d1 > 100);
                        jidx = find(jidx);
                        pidx = find((jidx - fidx) > 0);
                        nidx = find((jidx - fidx) < 0);
                        pdist = abs(jidx(pidx) - fidx);
                        ndist = abs(jidx(nidx) - fidx);
                        [~,spidx] = sort(pdist);
                        [~,snidx] = sort(ndist);
                        jidx = jidx([nidx(snidx(1)) pidx(spidx(1))]);
                    else
                        topJunctionPoint = obj.contourSequence{tm-1}.segment(:,obj.contourSequence{tm-1}.junctionIndex(1));
                        bottomJunctionPoint = obj.contourSequence{tm-1}.segment(:,obj.contourSequence{tm-1}.junctionIndex(2));
                        jidx(1) = obj.contourSequence{tm}.findNearestPoint(topJunctionPoint);
                        jidx(2) = obj.contourSequence{tm}.findNearestPoint(bottomJunctionPoint);
                    end
                    obj.contourSequence{tm}.tagJunctions(jidx);
                    %fprintf(['Done with tracing:' num2str(tm) ':' num2str(numel(obj.contourSequence)) '\n']);
                end
            
            catch ME
                ME;
            end
            
           
                
                
        end
        
        function [length] = getMidlineLength(obj)
            for tm = 1:numel(obj.contourSequence)
                length(tm) = obj.contourSequence{tm}.getMidlineLength();
            end
        end
        
        function [tipAngle] = getTipAngle(obj)
            for tm = 1:numel(obj.contourSequence)
                tipAngle(tm) = obj.contourSequence{tm}.getTipAngle();
            end
        end
        
        function [] = flipTraceDirectionToClockWise(obj)
            for e = 1:numel(obj.contourSequence)
                obj.contourSequence{e}.toClockWise();
            end
        end
    end
end