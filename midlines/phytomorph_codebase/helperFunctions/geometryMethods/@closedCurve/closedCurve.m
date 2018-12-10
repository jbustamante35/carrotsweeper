classdef closedCurve < matlab.mixin.Copyable & handle
    properties
        segment;            % curve
        filename;           % filename for the curve
        centerpoint;        % centerpoint for the curve
        
        skelVec;            % vector field pointing to the skeleton
        skelDist;           % distance to the skeleton
        
        centerVec;          % vector field pointing to the center
        centerDist;         % distance to the center point
        
        labels;             % labels
    end
    
    methods
        
        function [obj] = closedCurve(initObject)
            if isa(initObject,'closedCurve')
                obj = copy(initObject);
            else
                obj.segment = initObject;
                obj.centerpoint = mean(obj.segment,2);
                obj.calcDistaceToSkeleton();
                obj.calcDistanceToCenter();    
            end
        end
        
        function [] = setSegment(obj,segment)
            obj.segment = segment;
            obj.centerpoint = mean(obj.segment,2);
            obj.calcDistaceToSkeleton();
            obj.calcDistanceToCenter();
        end
        
        function [] = displace(obj,dx)
             newSegment = bsxfun(@plus,obj.segment,dx);
             obj.setSegment(newSegment);
             if ~isempty(obj.midline)
                 newMidline = bsxfun(@plus,obj.midline,dx);
                 obj.setMidline(newMidline);
             end
        end
        
        function [] = setMidline(obj,midline)
            obj.midline = midline;
        end
        
        function [] = calcDistaceToSkeleton(obj)
            [BW tmp] = obj.generateMask();
            tmp = tmp + 1;
            SKEL = bwmorph(BW,'skel',inf);            
            [bwd idx] = bwdist(SKEL);
            lidx = sub2ind(size(BW),tmp(2,:),tmp(1,:));            
            expD = bwd(lidx);
            quiv = idx(lidx);
            [q1 q2] = ind2sub(size(BW),quiv);
            Q = [q2;q1];
            Q = double(Q) - tmp;
            obj.skelVec = Q;
            obj.skelDist = expD;
        end
        
        function [] = calcDistanceToCenter(obj)
            obj.centerVec = bsxfun(@minus,obj.centerpoint,obj.segment);
            obj.centerDist = sum(obj.centerVec.*obj.centerVec,1).^.5;
        end
        
        function [sidx] = labelCurve(obj)
            [BW curve] = obj.generateMask();
            [r c V] = impixel(BW);
            obj.labels = zeros(numel(r),size(obj.segment,2));
            for e = 1:numel(r)
                dist = bsxfun(@minus,curve,[r(e);c(e)]);
                dist = sum(dist.*dist,1);
                [M sidx(e)] = min(dist);
                obj.labels(e,sidx(e)) = 1;
            end
            imshow(BW,[]);
            hold on
            plot(curve(1,:),curve(2,:),'r');
            plot(curve(1,sidx),curve(2,sidx),'g*');
            hold off;
            drawnow
        end
        
        function [BW tmp] = generateMask(obj)
            tmp = bsxfun(@minus,obj.segment,obj.centerpoint);
            mX = min(tmp,[],2);
            MX = max(tmp,[],2);
            LEN = round(MX - mX);
            tmp = bsxfun(@minus,tmp,mX);
            %tmp = tmp + 1;
            tmp = round(tmp);
            BW = poly2mask(tmp(1,:),tmp(2,:),LEN(2)+2,LEN(1)+2); 
        end
        
        function [I offSet] = getImage(obj,PAD)
            if nargin == 1;PAD = [0;0];end            
            offSet = obj.getOffSet - PAD;            
            BOX = obj.getBoundingBox();            
            ROWS = round([BOX(1)-PAD(1) BOX(1)+BOX(3)+PAD(1)]);
            COLS = round([BOX(2)-PAD(2) BOX(2)+BOX(4)+PAD(2)]);
            I = double(imread(obj.filename,'PixelRegion',{COLS ROWS}));
            offSet = [ROWS(1);COLS(1)];
        end
        
        function [BOX] = getBoundingBox(obj)
            UL = min(obj.segment,[],2);
            BR = max(obj.segment,[],2);
            BOX = [UL BR-UL]; 
        end
        
        function [BOX] = getPadBoundingBox(obj,PAD)
            BOX = obj.getBoundingBox();
            BOX(1:2) = BOX(1:2) - PAD;
            BOX(3:4) = BOX(3:4) + 2*PAD;
        end
        
        function [offSet] = getOffSet(obj)
            offSet = min(obj.segment,[],2);            
        end
        
        function [curveAtoms] = atomize(obj,NH,nhSZ,segmentSize)
            I = double(imread(obj.filename));
            curveAtoms = sampleCurveBank(I,obj,NH,nhSZ,segmentSize,@()closedCurveSegment());
            
            for e = 1:numel(curveAtoms)
                curveAtoms(e).centerVec = obj.centerpoint - curveAtoms(e).centerPoint;
                curveAtoms(e).centerDistance = dist(curveAtoms(e).centerVec);
                curveAtoms(e).skeletonVec = obj.skelVec(:,e);
                curveAtoms(e).skeletonDistance = obj.skelDist(e);
                curveAtoms(e).imageName = obj.filename;
                curveAtoms = obj.labelAtoms(curveAtoms);
            end
            
            
        end
        
        function [curveAtoms] = atomizeCurve(obj,segmentSize,containerClass)
            tmp = [obj.segment];
            for e1 = 1:size(tmp,2)
                % get curve segment
                sig = tmp(:,1:segmentSize)';
                % fit with least squares spline
                fn = spap2(1,3,[1:size(sig,1)]',sig');
                % get function derivative
                fn1 = fnder(fn,1);
                % eval the der for tangent vec
                tvec = fnval(fn1,(segmentSize-1)/2);
                tvec = tvec/norm(tvec);
                nvec = [tvec(2);-tvec(1)];
                % get the reference frame
                E = [tvec nvec];
                % get the center point
                U = fnval(fn,(segmentSize-1)/2);
                C = PCA_REPROJ(sig,E,U');            
                %{
                plot(U(1),U(2),'r.');
                quiver(U(1),U(2),E(1,2),E(2,2),'r');
                hold on
                drawnow
                %}
                tmp = circshift(tmp,[0 -1]);

                curveAtoms(e1) = containerClass();
                curveAtoms(e1).orientation = E;        
                curveAtoms(e1).segment = sig';
                curveAtoms(e1).Osegment = C';
                curveAtoms(e1).centerPoint = U;
            end
            curveAtoms = circshift(curveAtoms,[0 (segmentSize-1)/2]);
        end
        
        function [curveAtoms] = labelAtoms(obj,curveAtoms)
            if ~isempty(obj.labels)
                for e = 1:numel(curveAtoms)
                    curveAtoms(e).labelVec = obj.labels(:,e);
                end
            end
        end
        
        function [] = zeroCenterCurve(obj)
            obj.segment = bsxfun(@minus,obj.segment,obj.centerpoint);
        end
        
        function [] = displaceCurve(obj)
            obj.segment = bsxfun(@plus,obj.segment,obj.centerVec);
            obj.midline = bsxfun(@plus,obj.midline,obj.centerVec);
        end
        
        function [segments] = getIntersectionProfile(obj)
            N = closedCurve.getNormals(obj.segment,5);
            clearZone = 2;
            for e = 1:size(N,2)
                No = [obj.segment(:,e) obj.segment(:,e) - 2000*N(:,e)];
                SEG = [obj.segment];
                INT = InterX(SEG,No);
                tmp = bsxfun(@minus,INT,obj.segment(:,e));
                tmp = sum(tmp.*tmp,1).^.5;
                ridx = find(tmp < clearZone);
                if numel(ridx) == numel(tmp)
                    [J nridx] = max(tmp);
                    ridx(nridx) = [];
                end
                INT(:,ridx) = [];
                IS(:,e) = INT(:,end) - obj.segment(:,e);
                %{
                plot(obj.segment(1,:),obj.segment(2,:),'b');
                hold on
                quiver(obj.segment(1,e),obj.segment(2,e),IS(1,e),IS(2,e),0);
                hold off
                drawnow
                %}
            end
            segments = IS;
        end
        
        function [b] = contains(obj,target)
            ret = 0;
            b = inpoly(target.segment',obj.segment');
            b = all(b);
        end
        
        function [b] = containsPoint(obj,point)
            b = inpoly(point',obj.segment');
        end
        
        function [d] = hausdorfDistance(obj,target)
            for t = 1:numel(target)
                dist = [];
                for e = 1:size(obj.segment,2)
                    Q = bsxfun(@minus,target(t).segment,obj.segment(:,e));
                    dist(e) = min(sum(Q.*Q,1).^.5);
                end
                d(t) = mean(dist);
            end
        end
  
        function [SNIP] = getWindow(obj,windowLength,windowStart)
            SNIP = circshift(obj.segment,[0 -(windowStart-(windowLength-1)/2)]);
            SNIP = SNIP(:,1:(windowLength));
        end
        
        function [SNIP] = getAntiWindow(obj,windowLength,windowStart)
            SNIP = circshift(obj.segment,[0 -(windowStart-(windowLength-1)/2)]);
            SNIP = SNIP(:,windowLength+1:end);
        end
        
        function [] = flipTraceDirection(obj)
            tmp = flipdim(obj.segment,2);
            obj.setSegment(tmp);
        end
        
        function [] = toClockWise(obj)
            B = ispolycw(obj.segment(1,:),obj.segment(2,:));
            if ~B
                obj.flipTraceDirection();
            end
        end
    end
    
    methods (Static)
    
        function [K] = kurvature(segment,S)
            sz = size(segment);
            J = [segment';segment';segment'];

            % calculate curvature
            d1X1 = cwt(J(:,1),S,'gaus1');
            d1X2 = cwt(J(:,2),S,'gaus1');
            %d2X1 = cwt(J(:,1),S,'gaus2');
            %d2X2 = cwt(J(:,2),S,'gaus2');
            tan = cat(3,d1X1,d1X2);
            L = sum(tan.*tan,3).^.5;
            T = bsxfun(@times,tan,L.^-1);
            N = cat(3,T(:,:,2),-T(:,:,1));
            dT = gradient(T);
            K = sum(dT.*N,3);
            K = K.*L.^-1;
            %K = (d1X1.*d2X2 - d1X2.*d2X1).*(d1X1.^2 + d1X2.^2).^-3/2;            
            K = K(:,sz(2)+1:sz(2)+sz(2));
        end
    
        function [R] = distanceRatio(segment,S)
            
            CC = [segment segment segment];
            TAN = gradient(CC);
            L = sum(TAN.*TAN,1).^.5;
            L = im2col(L,[1 S],'sliding');
            L = sum(L,1);
            START = size(segment,2)-(S-1)/2+1;
            STOP  = START + size(segment,2) - 1;
            L = L(START:STOP);
            
            
            
            
            dCC = circshift(segment,-S);
            dCC = dCC - segment;
            dCC = sum(dCC.*dCC,1).^.5;
            
            R = L.*dCC.^-1;
        end    
    
        function [N] = getNormals(segment,S)
            sz = size(segment);
            J = [segment';segment';segment'];

            % calculate curvature
            d1X1 = cwt(J(:,1),S,'gaus1');
            d1X2 = cwt(J(:,2),S,'gaus1');            
            tan = cat(3,d1X1,d1X2);
            L = sum(tan.*tan,3).^.5;
            T = bsxfun(@times,tan,L.^-1);
            N = cat(3,T(:,:,2),-T(:,:,1));
            N = squeeze(N)';
            N = N(:,sz(2)+1:sz(2)+sz(2));
        end
        
        function [T] = getTangents(segment,S)
            sz = size(segment);
            J = [segment';segment';segment'];

            % calculate curvature
            d1X1 = cwt(J(:,1),S,'gaus1');
            d1X2 = cwt(J(:,2),S,'gaus1');            
            tan = cat(3,d1X1,d1X2);
            L = sum(tan.*tan,3).^.5;
            T = bsxfun(@times,tan,L.^-1);
            T = T(:,sz(2)+1:sz(2)+sz(2));
        end
        
        function [newC] = repara(segment)
            
        end
        
        function [w] = windowData(data,segmentSize)
            sz = size(data);
            data = [data data data];
            data = circshift(data,[0 (segmentSize-1)/2]);
            for e = 1:(size(data,2) - segmentSize)
                w(:,:,e) = data(:,e:e+(segmentSize-1));
            end
            w = w(:,:,1:sz(2));
        end
    end
        
end