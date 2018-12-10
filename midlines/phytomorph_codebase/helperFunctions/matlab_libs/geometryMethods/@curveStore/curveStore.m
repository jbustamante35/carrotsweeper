classdef curveStore < handle
    properties
        store;
    end
    
    methods
        
        function [obj] = curveStore()

        end
        
        function [] = insertParticles(obj,particles)
            for e = 1:numel(particles)
                obj.store{end+1} = particles(e);
            end
        end
        
    end
    
   
    methods (Static)
        
        function [b] = containment(source,targets)
            b = zeros(1,numel(targets));
            for e = 1:numel(targets)
                b(e) = source.contains(targets(e));
            end
        end
        
        function [b] = containmentMap(curves)
            b = zeros(size(curves));
            parfor e = 1:numel(curves)
                b(e,:) = curveStore.containment(curves(e),curves);
            end
            
            
            % find root nodes
            fidx = find(sum(b,1) == 1);
            % for each root node - it must contain objects
            for e = 1:numel(fidx)
                containmentValue(e) = sum(b(fidx(e),:),2);
            end
            
            
            
           %{
            for c = 1:numel(fidx)
                graph{c} = fidx(c);
                
                cidx = find(t(fidx(c),:));
                cidx = setdiff(cidx,fidx(c));
                d = curves(fidx(1)).hausdorfDistance(curves(cidx));
            end
            %}
        end
        
        function [children] = directChildren(t,nodeIndex)
            gc = [];
            children = curveStore.indirectChildren(t,nodeIndex);
            for c = 1:numel(children)
                gc = [gc curveStore.indirectChildren(t,children(c))];
            end
            children = setdiff(children,gc);
        end
        
        function [index] = indirectChildren(t,nodeIndex)
            index = setdiff(find(t(nodeIndex,:)),nodeIndex);
        end
        
        function [Q] = refine(innerCurve,outerCurve,image)
            newC = [];
            for e = 1:size(outerCurve.segment,2)
                delta = bsxfun(@minus,innerCurve.segment,outerCurve.segment(:,e));
                dis = sum(delta.*delta,1);
                [D(e) sidx] = min(dis);
                Q(:,e) = delta(:,sidx);
                X = linspace(outerCurve.segment(1,e),outerCurve.segment(1,e)+Q(1,e),100);
                Y = linspace(outerCurve.segment(2,e),outerCurve.segment(2,e)+Q(2,e),100);
                sample = ba_interp2(image,X,Y);
                [J midx] = max(sample);
                newC(:,e) = [X(midx);Y(midx)];
            end
        end
        
        function [fidx] = findRoots(containmentMap)
            fidx = find(sum(containmentMap,1) == 1);
        end
        
        function [curve] = fuseCurves(innerCurve,outerCurve,containmentClass)
            newC = [];
            for e = 1:size(outerCurve.segment,2)
                delta = bsxfun(@minus,innerCurve.segment,outerCurve.segment(:,e));
                dis = sum(delta.*delta,1);
                [D(e) sidx] = min(dis);                
                newC(:,e) = outerCurve.segment(:,e) + .5*delta(:,sidx);
            end            
            curve = containmentClass(newC);
            curve.filename = outerCurve.filename;
        end
        
        function [distance] = averageDistance(innerCurve,outerCurve)            
            if isempty(innerCurve) | isempty(outerCurve)
                distance = inf;
                return;
            end
            for e = 1:size(outerCurve.segment,2)
                delta = bsxfun(@minus,innerCurve.segment,outerCurve.segment(:,e));
                dis = sum(delta.*delta,1).^5
                [D(e) sidx] = min(dis);
            end
            distance = mean(D);
        end
        
        function [distance] = maxDistance(innerCurve,outerCurve)
            if isempty(innerCurve) | isempty(outerCurve)
                distance = inf;
                return;
            end
            for e = 1:size(outerCurve.segment,2)
                delta = bsxfun(@minus,innerCurve.segment,outerCurve.segment(:,e));
                dis = sum(delta.*delta,1).^.5;
                [D(e) sidx] = min(dis);
            end
            distance = max(D);
        end
        
        function [repCurve] = generateRepCurve(containmentMap,curves,threshold)
            fidx = curveStore.findRoots(containmentMap);
            for e = 1:numel(fidx)                
                parentIndex = fidx(e);
                % init rep curve
                repCurve(e) = curves(parentIndex);
                childIndex = curveStore.directChildren(containmentMap,parentIndex);
                if numel(childIndex) == 1
                    maxDistance = curveStore.maxDistance(curves(childIndex),curves(parentIndex));
                else
                    maxDistance = inf;
                end
                while numel(childIndex) == 1 & maxDistance < threshold
                    % fuse curves
                    repCurve(e) = curveStore.fuseCurves(curves(childIndex),repCurve(e),@(x)maizeContour(x));
                    %                     
                    parentIndex = childIndex;
                    childIndex = curveStore.directChildren(containmentMap,parentIndex);
                    if numel(childIndex) == 1
                        maxDistance = curveStore.maxDistance(curves(childIndex),curves(parentIndex));
                    else
                        maxDistance = inf;
                    end
                end
            end
        end
        
        function [curves] = containsPoints(curves,points,op)
            for e = 1:numel(curves)
                b = curves(e).containsPoint(points);
                b = op(b);
                if b
                    rm(e) = 0;
                else
                    rm(e) = 1;
                end
            end
            curves(find(rm)) = [];
        end
        
        function [] = deform(source,target)
            
            [sP sN] = source.toTensorForm();
            [tP tN] = target.toTensorForm();
            
            stepSize = .0001;         
            
            for loop = 1:10000
                delta = tP - sP;
                N = closedCurve.getNormals(sP,5);
                delta = sum(delta.*N,1);
                delta = bsxfun(@times,delta,N);
                delta = delta*stepSize;
                sP = sP + delta;
             
                plot(sP(1,900:1200),sP(2,900:1200),'r')
                hold on
                plot(tP(1,900:1200),tP(2,900:1200),'b')                
                drawnow      
                hold off
             
            end
                
        end
        
        function [delta] = potential(staticPoints,testPoint)
            delta = bsxfun(@minus,staticPoints,testPoint);
            distance = sum(delta.*delta,1).^.5;
            mag = distance.^-3;
            normalization = gamcdf(mag,1,2);
            mag = mag.*normalization;
            delta = bsxfun(@times,mag,delta);
            delta = sum(delta,2);
        end
        
        function [P] = curveProb(curveSet,curve)
            %{
            [vecForm normalSpace] = maizeContour.getExpectedCurve(curveSet);
            normalSpace = mean(normalSpace,3);            
            [iP iN] = curve.toTensorForm();
            ROT = sum(iN.*normalSpace,1);
            P = ROT;
            %}
            try
                
                
                % init vars
                downSample = 1;
                kScales = [5];
                tmpCurve = [];
                tmpNormal = [];
                K = [];
                % for each curve - get kurvature
                for e = 1:numel(curveSet)
                    [P N] = curveSet(e).toTensorForm();
                    
                    np = size(P,2);

                    P = P(:,1:downSample:end);
                    N = N(:,1:downSample:end);
                    K = cat(3,K,closedCurve.kurvature(P,kScales));
                    tmpCurve = cat(3,tmpCurve,P);
                    tmpNormal = cat(3,tmpNormal,N);
                end
                
                
                
                [iP iN] = curve.toTensorForm();
                fP = iP;
                iP = iP(:,1:downSample:end);
                iK = closedCurve.kurvature(iP,kScales);
                szK = size(iK);


                % compute prob
                K = reshape(K,[size(K,1)*size(K,2) size(K,3)]);
                iK = reshape(iK,[size(iK,1)*size(iK,2) size(iK,3)]);
                U = mean(K,2);
                S = std(K,1,2);
                %for e = 1:size(U,1)
                %    f(e) = normpdf(iK(e),U(e),S(e));
                %end
                Z = ((iK(:,:) - U(:,:)).*S(:,:).^-1);
                Z = reshape(Z,szK);
                f = normpdf(iK,U,S);
                %f = f.*mean(S(:));
                %max = normpdf(U,U,S);
                
                %f = f.*max.^-1;
                %{
                for e = 1:size(K,1)
                    [f(e) xi(e) u(e)] = ksdensity(K(e,:),iK(e));
                    f(e) = f(e)*u(e);
                end
                %}
                
                
                
                f = reshape(f,szK);
                S = reshape(S,szK);
                uS = mean(S,2);
                f = bsxfun(@times,f,uS);
                
                f = mean(f,1);
                %f = interp1(1:numel(f),f,linspace(1,numel(f),np));
                P = f;
                
                thresh = .05;
                if any(P < thresh)
                    here = 1;
                    fidx = find(P < thresh);
                    plot(fP(1,:),fP(2,:),'r');
                    hold on
                    plot(fP(1,fidx),fP(2,fidx),'k*');
                    axis equal
                    waitforbuttonpress;
                end
                
                R = regionprops(P < thresh,'pixelIdx');
                
                % generate curve patch
                
                
                
            catch ME
                P = 0;
            end
            
        end
        
        function [P] = createCurvePatch(curveSet,window,point)                                    
            tmpCurve = [];
            tmpNormal = [];
            for e = 1:numel(curveSet)
                [P N] = curveSet(e).toTensorForm();                
                tmpCurve = cat(3,tmpCurve,P);
                tmpNormal = cat(3,tmpNormal,N);
            end
        end
    end
end