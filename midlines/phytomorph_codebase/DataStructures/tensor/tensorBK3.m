classdef tensor < handle
    properties
        M;
    end
    
    methods
        % constructor
        function [obj] = tensor(M)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting construction operation:\n']);tm = clock;
            % default assign the tensor to []
            if nargin == 0
                M = [];
            end
            % if the first object isa tensor -  then init properly
            if isa(M,'tensor')
                obj.M = M.M;
            else
                % remove the singleton dims
                M = squeeze(M);
                % transpose/permute - rows and column vectors are not things in this universe
                if size(M,1) == 1
                    M = permute(M,[2 1]);
                end
                % convert to double if needed
                if ~isa(M,'double')
                    M = double(M);
                end
                % assign
                obj.M = M;
            end
            fprintf(['\tending construction operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
        
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'.')
                try 
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                catch
                    try
                        [varargout{1:nargout}] = builtin(S(1).subs,obj.M,S(2).subs{:});
                        for e = 1:numel(nargout)
                            varargout{e} = tensor(varargout{e});
                        end
                    catch
                         fprintf(['error in index-. for tensor.\n']);
                    end
                end
            elseif strcmp(S(1).type,'()')
                [varargout{1:nargout}] = builtin('subsref',obj.M,S);
            end
        end
        
        % get the order of the tensor
        function [r] = order(obj)
            r = size(obj.M);
            r(r==1) = [];
            r = numel(r);
        end
        
        
        % fold the tensor object
        function [] = fold(obj,foldDirections)
            % decode the fold instructions into a permute instruction operand
            permuteOperand = [];
            oldSZ = size(obj.M);
            newSZ = zeros(1,numel(foldDirections));
            for e = 1:numel(foldDirections)
                % this builds the permute operand
                permuteOperand = [permuteOperand foldDirections{e}];
                % post permute this will contain the needed size information
                newSZ(e) = [prod(oldSZ(foldDirections{e}))];
            end
            % if the newSZ is a column/row vector then make [column 1]
            if numel(newSZ) == 1
                newSZ = [newSZ 1];
            end
            if numel(permuteOperand) == 1
                permuteOperand = [permuteOperand 2];
            end
            % check that all dims are accounted for
            dF = setdiff(min(permuteOperand):max(permuteOperand),1:numel(oldSZ));
            if isempty(dF)
                obj.M =  permute(obj.M,permuteOperand);
                obj.M = reshape(obj.M,newSZ);
            else
                fprintf(['error in fold operation due to lack of dims accounting.\n']);
            end
        end
        
        
        function [r] = or(objA,objB)
            if objB.order == 2 && objA.order == 1
                r = tensor.dotProduct(objA,objB,[1 1]);
            elseif objA.order == 2 && objB.order == 1
                r = tensor.dotProduct(objA,objB,[2 1]);
            else
                r = tensor.dotProduct(objA,objB,[1 1]);
            end
           
        end
        
        function [r] = gt(objA,objB)
            objA.fold({[1:objA.order()]});
            r = objA;
        end
        
        function [r] = lt(objA,objB)
            objB.fold({[1:objB.order()]});
            r = objB;
        end
    end
    
    
    methods (Static)
        
        
        
         % apply dot without respect the the typing of the tensor
        function [r] = dotProduct(a,b,d)
            if tensor.dimsCheck(a.M,b.M,d)
                [a,sza,newOrder] = tensor.pullFront(a.M,d(1));
                [b,szb,newOrder] = tensor.pullFront(b.M,d(2));
                r = mtimesx(a,'T',b);
                r = reshape(r,[sza(2:end) szb(2:end)]);
                r = tensor(r);
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
            if nargin == 2
                if tensorContainer.orderCheck(a,b)
                    d = [1:a.order() ; 1:b.order()]';
                else
                    fprintf(['order does not match therefore dims can not!\n']);
                    c = false;
                end
            end
            szA = size(a);
            szB = size(b);
            
           
            ca = szA(d(:,1));
            cb = szB(d(:,2));
            c = all(all(ca == cb,2),1);
        end
        
        function [c] = orderCheck(a,b)
            c = a.order() == b.order();
        end
        
        
        function [array,sz,newOrder] = pullFront(array,order)
            newOrder = [order setdiff(1:ndims(array),order)];
            array = permute(array,newOrder);
            sz = size(array);
            array = reshape(array,[sz(1) prod(sz(2:end))]);
        end
        
        function [array] = pushBack(array,sz,order)
            array = reshape(array,sz);
            array = ipermute(array,order);
        end
        
    end
end






