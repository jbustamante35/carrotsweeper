classdef phytoPath < handle
    properties
        pathData;
        childList;
        parent;
        RSAblock;
    end
    
    methods
        
        function [obj] = phytoPath(RSA,varargin)
            obj.RSAblock = RSA;            
            if nargin == 2
                obj.pathData = varargin{1};
            end
        end
        
        function [] = setPathData(obj,pathData)
            obj.pathData = pathData;
        end
        
        function [] = attachChild(obj,child,location)
            obj.childList(end+1).child = child;
            obj.childList(end).location = location;
            obj.RSAblock.addBranch(child);
            child.parent = obj;
        end
        
        function [] = plot(obj,varargin)
            L = 0;
            NM = 0;
            CL = 'b';  
            if nargin == 2
                CL = varargin{1};                
            elseif nargin >2
                L = varargin{2};
                NM = varargin{3};            
            end
            SNIP = 3;
            [A E] = measureTipAngle(obj,SNIP);
            plot(obj.pathData(2,:),obj.pathData(1,:),CL);
            %textLocation = mean(obj.pathData,2);
            textLocation = mean(obj.pathData(:,end),2);
            if L ~= 0 & NM ~= 0                
                text(textLocation(2),textLocation(1),[num2str(L) '.' num2str(NM)],'Background','w','FontUnits','pixels','FontSize',1);
            end
            quiver(obj.pathData(2,end),obj.pathData(1,end),E(2),E(1),10,'Color','r');
            
            for ch = 1:numel(obj.childList)                
                obj.childList(ch).child.plot('g',L+1,ch);
            end
        end
        
        function [dL] = measureLength(obj)
            dL = diff(obj.pathData,1,2);
            dL = sum(sum(dL.*dL,1).^.5);
        end
        
        function [A] = measureSkew(obj)
            vec = obj.pathData(:,end) - obj.pathData(:,1);
            A = atan2(vec(2),vec(1));
        end
        
        function [A E] = measureTipAngle(obj,SNIP)
            % tip angle
            str = max(size(obj.pathData,1)-SNIP+1,1);
            % pca on SNIP
            [S C U E L ERR LAM] = PCA_FIT_FULL(obj.pathData(:,str:end)',1);
            % orient via diff
            pD = diff(obj.pathData(:,str:end),1,2);
            pD = mean(pD,2);
            if pD'*E < 0
                E = -E;
            end        
            % return values
            E = E;
            A = acos(E(1)/norm(E))*180/pi;
        end
        
        function [K] = measureKurvature(obj,SM)
            K = cwtK_filter(obj.pathData',{SM});
        end
        
        function [H1 H2] = measurePhenotypes(obj,varargin)
            if nargin == 1
                H1 = {'Branch' 'Length' 'Skew' 'Tip Angle'};
                LV = 0;
                BN = 0;
                H2 = {};
            else
                H1 = varargin{1};
                LV = varargin{3};
                BN = varargin{4};
                H2 = varargin{2};
            end
            
            L = obj.measureLength();
            SW = obj.measureSkew();
            TA = obj.measureTipAngle(20);
            B = [num2str(LV) '.' num2str(BN)];
            
            H1{end+1,1} = B;
            H1{end,2} = L;
            H1{end,3} = SW;
            H1{end,4} = TA;
            
            H2{1,end+1} = B;
            K = obj.measureKurvature(5);
            for e = 1:numel(K.K)
                H2{e+1,end} = K.K(e);
            end
            
            for e = 1:numel(obj.childList)
                [H1 H2] = obj.childList(e).child.measurePhenotypes(H1,H2,LV+1,e);
               
            end
            
        end
    end
end