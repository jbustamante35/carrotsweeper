classdef grt <  matlab.mixin.Heterogeneous
    % general recursive tensor
    properties
        T;
    end
    
    
    methods
    
    end
    
    methods (Static)
        
        function [r] =  map(objA,objB,d,op)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % while the map operation (op) is defined along the 
            % dim array d, the operation takes each element and
            % therefore all dims must be equal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % note that if objects are symbols then need to call different constructor
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if operands are numeric - then apply directly
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isnumeric(objA) && isnumeric(objB)
                r = op(objA,objB);
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if operands are not arrays of many - then deep dive
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if numel(objA) == 1 && numel(objB) == 1
                    % pop/grab dim instructions from d
                    newDimI = d(2:end,:);
                    if grt.dimsCheck(objA,objB)
                        for e = 1:numel(objA)
                            r(e) = tensor(grt.map(objA(e).T,objB(e).T,op,newDimI));
                        end
                        r = reshape(r,size(objA));
                    end
                else
                    % check if objects have the same dim
                    % resuse the dim instructions until deep dive
                    if grt.dimsCheck(objA,objB)
                        for e = 1:numel(objA)
                            r(e) = grt.apply(objA(e),objB(e),op,d);
                        end
                    end
                    r = reshape(r,size(objA));
                    [r] = grt.fixOrder(r);
                end
            end
            
            
            
            if grt.dimsCheck(objA,objA)
                [objA,szA,newOrderA] = pullFront(objA,d(1));
                [objA,szB,newOrderB] = pullFront(objB,d(2));
                % due to dims check the front dims are now equal
                for e1 = 1:size(A,1)
                    for e2 = 1:size(A,2)
                        r(e1,e2) = op(objA(e1,e2),objB(e1,e2));
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % recursivly apply operation op with objA and objB as the operands
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [r] = apply(objA,objB,op)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if operands are numeric - then apply directly
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if isnumeric(objA) && isnumeric(objB)
                r = op(objA,objB);
            else
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if operands are not arrays of many - then deep dive
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if numel(objA) == 1 && numel(objB) == 1
                    if grt.dimsCheck(objA,objB)
                        for e = 1:numel(objA)
                            r(e) = tensor(grt.apply(objA(e).T,objB(e).T,op));
                        end
                        r = reshape(r,size(objA));
                    end
                else
                    % same
                    if grt.dimsCheck(objA,objB)
                        for e = 1:numel(objA)
                            r(e) = grt.apply(objA(e),objB(e),op);
                        end
                    end
                    r = reshape(r,size(objA));
                    [r] = grt.fixOrder(r);
                end
            end
           
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % apply dot without respect the the typing of the tensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [r] = dotProduct(a,b,d)
            
            
            if grt.dimsCheck(a,b,d(1,:))
                
                [a,sza,newOrder] = grt.pullFront(a,d(1));
                [b,szb,newOrder] = grt.pullFront(b,d(2));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if the operands are numeric - need to handle symbolic too
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if isnumeric(a) && isnumeric(b)
                    r = mtimesx(a,'T',b);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % else if the operands are objects that are holding numeric/symbolic
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % assume that size(a,1) == size(a,2) 
                    % should be by defined above checks
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for e = 1:size(a,1)
                        r(e) = grt.dotProduct(a(e).T,b(e).T,d(2:end,:));
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % reshape to the proper size
                r = reshape(r,[sza(2:end) szb(2:end)]);
                r = tensor(r);
            else
                fprintf(['dotProduct not allowed\n']);
            end
            
            
            
        end
        
        
        function [r] = tensorProduct(a,b)
            r = mtimesx(a.T(:),b.T(:),'T');
            r = reshape(r,[size(a.T) size(b.T)]);
            r = tensor(r);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % general dims check 
        % if d is given then check those dims
        % else all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [c] = dimsCheck(a,b,d)
            %%%%%%%%%%%%%%%%%%%
            % if dims are not given - then dims check all
            %%%%%%%%%%%%%%%%%%%
            if nargin == 2
                if grt.orderCheck(a,b)
                    d = [1:grt.order(a) ; 1:grt.order(b)]';
                else
                    fprintf(['order does not match therefore dims can not!\n']);
                    c = false;
                end
            end
            %%%%%%%%%%%%%%%%%%%
            % get the size of the tensors
            szA = size(a);
            szB = size(b);
            
           
            ca = szA(d(:,1));
            cb = szB(d(:,2));
            % check
            c = all(all(ca == cb,2),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generalize order check at the level of query
        % will NOT yet recurse levels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [c] = orderCheck(a,b)
            c = grt.order(a) == grt.order(b);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the order of the tensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [r] = order(obj)
            tmp = size(obj);
            tmp(tmp==1) = [];
            r = numel(tmp);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pull to the "order" to the front and streak the rest
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [array,sz,newOrder] = pullFront(array,order)
            % get the new order of the result tensor
            newOrder = [order setdiff(1:ndims(array),order)];
            % call permute
            array = permute(array,newOrder);
            % get the size
            sz = size(array);
            % reshape with front face as order
            array = reshape(array,[sz(1) prod(sz(2:end))]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inverse of pullFront
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [array] = pushBack(array,sz,order)
            array = reshape(array,sz);
            array = ipermute(array,order);
        end
        
        function [T] = fixOrder(T)
            % remove the singleton dims
            T = squeeze(T);
            % transpose/permute - rows and column vectors are not things in this universe
            if size(T,1) == 1
                T = permute(T,[2 1]);
            end
        end
        
        
        function [r] = generalApply(varargin)
            inputs = varargin{1:end-2};
            % second last spot - function
            % last spot - arguments
            if isa(varargin{end-1},'sym')
                body = varargin{end-1};
                args = varargin{end};
                body(args) = body;
                
            elseif isa(varargin{end-1},'char')
                
            elseif isa(varargin{end-1},'function_handle')
                body = varargin{end-1};
                
            end
            r = body(inputs);
        end
        
        
        %{
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
        %}
        
        
        
        
        
        
        
        
        
        
    end
end