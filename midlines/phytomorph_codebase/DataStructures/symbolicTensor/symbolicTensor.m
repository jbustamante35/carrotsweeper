classdef symbolicTensor < grt
    
    properties
        order;
        c_interface;
        %forwardMap;
        %inverseMap;
    end
    
    methods
        function [obj] = symbolicTensor(map,order)
            obj.T = formula(map);
            obj.c_interface = argnames(map);
            %{
            if nargin == 2
                inverseMap = [];
            end
            %}
            
            %obj.forwardMap = map;
            %obj.inverseMap = inverseMap;
            obj.order = order;
        end
        
        
        function [] = apply(obj,varargin)
            varargin{end+1} = 'placeHolder';
            varargin{end+1} = 'placeHolder';
            for e = 1:numel(obj)
                varargin{end-1} = obj.T;
                varargin{end} = obj.c_interface;
                r = grt.generalApply(varargin{:})
            end
        end
        
        %{
        function [] = tryStuff(obj)
             af = formula(obj.forwardMap);
             sz = size(af);
             newOrder = [2 3 1 4];
             af = permute(af,newOrder);
             szNew = size(af);
             af = reshape(af,[szNew(1:2) prod(szNew(3:end))]);
             for e = 1:size(af,3)
                af2(:,:,e) = inv(af(:,:,e));
             end
             af2 = reshape(af2,szNew);
             af2 = ipermute(af2,newOrder);
        end
        %}
        
        %{
        function [r] = dot(objA,objB,d)
            [r] = symbolicTensor.dotProduct(objA,objB,d);
        end
        %}
        
        
        %{
        function [varargout] = applyForward(obj,varargin)
            
            
            toApply = obj(1);
            for e = 2:numel(obj)
                toApply = toApply.dot(obj(e),[1 1]);
            end
            
            
            
            outClass = str2func(class(varargin{1}));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % gather inputs
            
            nm = argnames(toApply.forwardMap);
            formulaMap = formula(toApply.forwardMap);
            cnt = 1;
            for inputNumber = 1:numel(varargin)
                toCatAlong(inputNumber) = ndims(varargin{inputNumber}.M);
                sz{inputNumber} = size(varargin{inputNumber}.M);
                tmpy = reshape(varargin{inputNumber}.M,[prod(sz{inputNumber}(1:end-1)) sz{inputNumber}(end)]);
                for value = 1:size(tmpy,2)
                    V{cnt} = tmpy(:,value);
                    cnt = cnt + 1;
                end
                nSZ(inputNumber) = prod(sz{inputNumber}(1:end-1));
            end
            [~,midx] = max(nSZ);
            outSize = sz{midx}(1:end-1);
            toCatAlong = max(toCatAlong);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % apply inputs 
            r = [];
            for outputNumber = 1:size(formulaMap,1)
                tmpFormula = matlabFunction(formulaMap(outputNumber),'Vars',nm);
                tmpOut = tmpFormula(V{:});
                tmpOut = reshape(tmpOut,outSize);
                r = cat(toCatAlong,r,tmpOut);
            end

            varargout{1} = outClass('','',r);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        
        function [varargout] = applyInverse(obj,varargin)
            
            obj = flipdim(obj,2);
            toApply = obj(1);
            for e = 2:numel(obj)
                toApply = toApply.dot(obj(e),[1 1]);
            end
            outClass = str2func(class(varargin{1}));

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % gather inputs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nm = argnames(toApply.inverseMap);
            formulaMap = formula(toApply.inverseMap);
            cnt = 1;
            for inputNumber = 1:numel(varargin)
                toCatAlong(inputNumber) = ndims(varargin{inputNumber}.M);
                sz{inputNumber} = size(varargin{inputNumber}.M);
                tmpy = reshape(varargin{inputNumber}.M,[prod(sz{inputNumber}(1:end-1)) sz{inputNumber}(end)]);
                for value = 1:size(tmpy,2)
                    V{cnt} = tmpy(:,value);
                    cnt = cnt + 1;
                end
                nSZ(inputNumber) = prod(sz{inputNumber}(1:end-1));
            end
            [~,midx] = max(nSZ);
            outSize = sz{midx}(1:end-1);
            toCatAlong = max(toCatAlong);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % apply inputs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            r = [];
            for outputNumber = 1:size(formulaMap,1)
                tmpFormula = matlabFunction(formulaMap(outputNumber),'Vars',nm);
                tmpOut = tmpFormula(V{:});
                tmpOut = reshape(tmpOut,outSize);
                r = cat(toCatAlong,r,tmpOut);
            end
            varargout{1} = outClass('','',r);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %}
    end
    
    
    
    methods (Static)
        %{
        % apply dot without respect the the typing of the tensor
        function [r] = dotProduct(a,b,d)
            if symbolicTensor.dimsCheck(a,b,d)
                %%%%%%%%%%%%%%%%%%%%%
                % make forward map
                argA = argnames(a.forwardMap);
                argB = argnames(b.forwardMap);
                af = formula(a.forwardMap);
                bf = formula(b.forwardMap);
                [af,sza,newOrderA] = tensor.pullFront(af,d(1));
                [bf,szb,newOrderB] = tensor.pullFront(bf,d(2));
                r = af*bf; 
                r = reshape(r,[sza(2:end) szb(2:end)]);
                newOrder = size(r);
                r([argB(:);argA(:)]) = r;
                r = symbolicTensor(r,newOrder,[]);
            else
                fprintf(['dotProduct not allowed\n']);
            end
        end
        
        function [r] = tensorProduct(a,b)
            r = mtimesx(a.M(:),b.M(:),'T');
            r = reshape(r,[size(a.M) size(b.M)]);
            r = tensor(r);
        end
        
        function [c] = dimsCheck(a,b,d)
           c = true;
        end
        %}
    end
    
    
    
    
    
    
    
    
    
end