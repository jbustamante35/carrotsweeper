classdef tensor < grt
    properties
        
    end
    
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = tensor(T)
            fprintf(['----------------------------------------------------------\n']);
            fprintf(['\tstarting construction operation:\n']);tm = clock;
            % default assign the tensor to []
            if nargin == 0;T = [];end
            % call to fix order
            [T] = grt.fixOrder(T);
            % if the first object isa tensor -  then init properly
            if ~isa(T,'tensor');T = double(T);end
            obj.T = T;
            fprintf(['\tending construction operation:' num2str(etime(clock,tm)) ':\n']);
            fprintf(['----------------------------------------------------------\n\n']);
        end
        
        %{
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
        %}
    end
    
    %{
    methods (Static)
        
        
        
        
        % fold the tensor object
        function [r] = apply(objA,objB,op)
            
            if isnumeric(objA) && isnumeric(objB)
                r = op(objA,objB);
            else
                % if given one object - then deep dive
                if numel(objA) == 1 && numel(objB) == 1
                    if tensor.dimsCheck(objA,objB)
                        for e = 1:numel(objA)
                            r(e) = tensor(tensor.apply(objA(e).M,objB(e).M,op));
                        end
                        r = reshape(r,size(objA));
                    end
                else
                    % same
                    if tensor.dimsCheck(objA,objB)
                        for e = 1:numel(objA)
                            r(e) = tensor.apply(objA(e),objB(e),op);
                        end
                    end
                    r = reshape(r,size(objA));
                    [r] = tensor.fixOrder(r);
                end
            end
           
            
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % apply dot without respect the the typing of the tensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [r] = dotProduct(a,b,d)
            
            
            if tensor.dimsCheck(a,b,d(1,:))
                
                [a,sza,newOrder] = tensor.pullFront(a,d(1));
                [b,szb,newOrder] = tensor.pullFront(b,d(2));
                
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
                        r(e) = tensor.dotProduct(a(e).M,b(e).M,d(2:end,:));
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
            r = mtimesx(a.M(:),b.M(:),'T');
            r = reshape(r,[size(a.M) size(b.M)]);
            r = tensor(r);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % general dims check 
        % if d is given then check those dims
        % else all
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [c] = dimsCheck(a,b,d)
            % if 
            if nargin == 2
                if tensor.orderCheck(a,b)
                    d = [1:tensor.order(a) ; 1:tensor.order(b)]';
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
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generalize order check at the level of query
        % will NOT yet recurse levels
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [c] = orderCheck(a,b)
            c = tensor.order(a) == tensor.order(b);
        end
        
        % get the order of the tensor
        function [r] = order(obj)
            tmp = size(obj);
            tmp(tmp==1) = [];
            r = numel(tmp);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pull to the "order" to the front and streak the rest
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [array,sz,newOrder] = pullFront(array,order)
            % 
            newOrder = [order setdiff(1:ndims(array),order)];
            array = permute(array,newOrder);
            sz = size(array);
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
    end
    %}
end






