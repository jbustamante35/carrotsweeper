classdef maizeContour < closedCurve
    properties
        totalGradient;
        tip;
        tipIndex;
        
        junction;
        junctionIndex;
        
        junctionMidpoint;
        
        midline;
        midlineFrame;
    end
    
    methods
        
        function [obj] = maizeContour(initObject)
            obj = obj@closedCurve(initObject);
        end
        
        function [] = tracemidline(obj)
            maxLength = 1500;
            midlineThresh = 10;
            
            [BW tmp] = obj.generateMask();
            sDIST = bwdist(~BW);             
            [I offSet] = obj.getImage([20;20]);
            I = imfilter(I,fspecial('gaussian',[11 11],5));
            [g1 g2] = gradient(I);
            newTip = obj.tip - offSet;
            t1 = ba_interp2(g1,newTip(1),newTip(2));
            t2 = ba_interp2(g2,newTip(1),newTip(2));
            N = [t2 t1];
            N = -N / norm(N);
            T = [N(2) -N(1)];
            initD = [T;N];
            tmpTip = obj.tip - obj.getOffSet();
            obj.midline = trackFromPointAtGradient(double(sDIST),tmpTip,initD,maxLength,15,pi,[15 200],.3);
            %distAlongMidline = ba_interp2(double(sDIST),obj.midline(1,:),obj.midline(2,:));
            %{
            fidx = distAlongMidline < midlineThresh;
            gidx = find(fidx==0);
            fidx = 1:gidx(1)-1;
            obj.midline = obj.midline(:,fidx);
            %}
            obj.midline = bsxfun(@plus,obj.midline,obj.getOffSet());
        end
        
        function [] = tagTip(obj,tipIndex)
            obj.tipIndex = tipIndex;
            obj.tip = obj.segment(:,obj.tipIndex);
            obj.tracemidline();
        end
     
        function [sidx] = labelCurve(obj)
            sidx = labelCurve@closedCurve(obj);
            obj.tagTip(sidx(1));
            obj.tagJunctions(sidx(2:3));
        end
        
        function [] = tagJunctions(obj,junctionIndex)
            obj.junctionIndex = junctionIndex;
            obj.junction = obj.segment(:,obj.junctionIndex);
            
            TAN = diff(obj.junction,1,2);
            TAN = TAN/norm(TAN);
            U = mean(obj.junction,2);
            NOR = [TAN(2) -TAN(1)]';
            E = [TAN NOR];
            comp = PCA_REPROJ(obj.midline',E,U');
            fidx = find(comp(:,2) > 0);
            
            obj.midline = obj.midline(:,1:fidx(1));
        end
        
        function [NOR] = generateNormalField(obj)
            NOR = [];
            segmentSize = 31;
            tmp = [obj.segment];
            for e1 = 1:size(tmp,2)
                % get curve segment
                sig = tmp(:,1:segmentSize)';
                % fit with least squares spline
                fn = spap2(1,2,[1:size(sig,1)]',sig');
                % get function derivative
                fn2 = fnder(fn,1);
                pvec = fnval(fn,(segmentSize-1)/2);
                % eval the der for tangent vec
                tvec = fnval(fn2,(segmentSize-1)/2);
                tvec = tvec/norm(tvec);
                nvec = [tvec(2);-tvec(1)];
                tmp = circshift(tmp,[0 -1]);
                NOR = [NOR -nvec];
            end
            NOR = circshift(NOR,[0 (segmentSize-1)/2]);
        end
       
        function [vecForm normalSpace kernel root] = toTensorForm(obj)
            RootSample = 501;
            kernelDownSample = 10;
            X = obj.junctionIndex(2):obj.tipIndex;
            Xi = linspace(X(1),X(end),501);
            Y = obj.segment(:,X);
            topRoot = interp1(X',Y',Xi')';            
            topRoot(:,end) = [];
            
            X = obj.tipIndex:obj.junctionIndex(1);
            Xi = linspace(X(1),X(end),500);
            Y = obj.segment(:,X);
            bottomRoot = interp1(X',Y',Xi')';
            
            root = [topRoot bottomRoot];
            
            tmp = obj.segment(:,1:end-1);
            tmp = circshift(tmp,[0 -obj.junctionIndex(1)]);
            snip = size(tmp,2) - obj.junctionIndex(1);
            X = 1:(obj.junctionIndex(2) + snip);
            Xi = linspace(X(1),X(end),1001);
            Y = tmp(:,X);
            kernel = interp1(X',Y',Xi')';
            kernel(:,end) = [];
            
            h = size(kernel,2)/2;
            half1 = kernel(:,h:end);
            half2 = kernel(:,1:h-1);
            
            vecForm = [half1 topRoot bottomRoot half2];
            
            %T = diff(vecForm,1,2);
            %L = sum(T.*T,1).^.5;
            %vecForm = interp1([0 cumsum(L)]',vecForm',linspace(L(1),sum(L),size(vecForm,2)));
            %vecForm = vecForm';
            
            [S C U E L ERR LAM] = PCA_FIT_FULL(vecForm',2);
            if sign(E(:,1)'*vecForm(:,obj.tipIndex)) == -1
                E(:,1) = -E(:,1);
            end
            E(:,2) = [E(2,1) -E(1,1)];
            
            vecForm = PCA_REPROJ(vecForm',E,U);
            vecForm = vecForm';
            
            [T] = gradient(vecForm);
            L = sum(T.*T,1).^-.5;
            T = bsxfun(@times,T,L);
            normalSpace = [T(2,:);-T(1,:)];
            kernel = kernel(:,1:kernelDownSample:end);
        end
        
        function [vecForm sz] = toVectorForm(obj)
            [position normal] = obj.toTensorForm();
            
            sz.position = [size(position)];
            sz.normal = [size(normal)];
            vecForm = [position(:)' normal(:)'];
            sz.bundle = size(vecForm);
        end
        
        function [curve] = centeredForm(obj)
            
        end
        
        function [kernelCurve] = getKernelCurve(obj)
            [tensorForm] = obj.toTensorForm();
        end
        
        function [] = plot(obj)
            lineWidth = 1;
            plot(obj.segment(1,:),obj.segment(2,:),'b','LineWidth',lineWidth);
            if ~isempty(obj.midline)
                hold on
                plot(obj.midline(1,:),obj.midline(2,:),'r','LineWidth',lineWidth);
            end
            if ~isempty(obj.junctionIndex)
                plot(obj.segment(1,obj.junctionIndex(1)),obj.segment(2,obj.junctionIndex(1)),'g*');
                plot(obj.segment(1,obj.junctionIndex(2)),obj.segment(2,obj.junctionIndex(2)),'g*');
                plot(obj.segment(1,obj.junctionIndex),obj.segment(2,obj.junctionIndex),'g','LineWidth',lineWidth);
            end
             if ~isempty(obj.tipIndex)
                plot(obj.segment(1,obj.tipIndex(1)),obj.segment(2,obj.tipIndex(1)),'k*');
             end
             hold off
        end
        
        function [distance] = meanDistance(obj,target)
            tmp1 = bsxfun(@minus,obj.segment,obj.centerpoint);
            tmp2 = bsxfun(@minus,target.segment,target.centerpoint);
            for pt = 1:size(obj.segment,2)
                delta = bsxfun(@minus,tmp1(:,pt),tmp2);
                delta = sum(delta.*delta,1).^.5;
                distance(pt) = min(delta);
            end
            distance = mean(distance);
        end
        
        function [fidx] = findNearestPoint(obj,point)
            delta = bsxfun(@minus,obj.segment,point);
            delta = sum(delta.*delta,1).^.5;
            [J,fidx] = min(delta);
        end
        
        function [length] = getMidlineLength(obj)
            if ~isempty(obj.midline)
                dl = diff(obj.midline,1,2);
                dl = sum(dl.*dl,1).^.5;
                length = sum(dl);
            end
        end
        
        function [tipAngle] = getTipAngle(obj)
            data = [obj.segment obj.midline];
            data = bsxfun(@minus,data,obj.tip);
            dL = sum(data.*data,1).^.5;
            fidx = find(dL < 25);
            [S C U E L ERR LAM] = PCA_FIT_FULL(data(:,fidx)',2);
            %{
            plot(data(1,:),data(2,:));
            hold on
            quiver(0,0,E(1,1),E(2,1),100,'g');
            quiver(0,0,N(1,obj.tipIndex),N(2,obj.tipIndex),100,'g');
            %}
            N = closedCurve.getNormals(obj.segment,5);
            if N(:,obj.tipIndex)'*E(:,1) < 0
                E(:,1) = -E(:,1);
            end
            %tipAngle = atan2(N(2,obj.tipIndex),N(1,obj.tipIndex));
            tipAngle = atan2(E(2,1),E(1,1));
        end
        
        
        
    end
    
    methods (Static)
        
        
        function [vecForm normalSpace] = getExpectedCurve(curveSet)
            for e = 1:numel(curveSet)
                [vecForm(:,:,e) normalSpace(:,:,e)] = curveSet(e).toTensorForm();
            end
        end
        
        function [] = plotCurveSet(curveSet)
            [position normal] = maizeContour.getExpectedCurve(curveSet);
            P = mean(position,3);
            N = mean(normal,3);
            quiver(P(1,:),P(2,:),N(1,:),N(2,:))
        end
        
        function [B] = decomposeCurveSet(curveSet)
            for e = 1:numel(curveSet)
                [bundle(:,e) sz] = curveSet(e).toVectorForm();
            end
            pSTART = 1;
            pSTOP = prod(sz.position);
            nSTART = pSTOP+1;
            nSTOP = nSTART + prod(sz.normal) - 1;
            [B.pS B.pC B.pU B.pE B.pL B.pERR B.pLAM] = PCA_FIT_FULL(bundle(pSTART:pSTOP),4);
            [B.nS B.nC B.nU B.nE B.nL B.nERR B.nLAM] = PCA_FIT_FULL(bundle(nSTART:nSTOP),4);
        end
    end
end
