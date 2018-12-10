classdef tensorContainer < matlab.mixin.Heterogeneous
    properties
        ten;
    end
    
    methods
        function [obj] = tensorContainer(ten)
            if nargin == 0
                ten = tensor(0);
            end
            if isa(ten,'tensorContainer')
                obj = ten;
            else
                obj.ten = ten;
            end
            
        end
        
        function [r] = plus(objA,objB)
            cnt = 1;
            % perform size check
            if all(size(objA) == size(objB))
                r(numel(objA)) = tensorContainer();
                for e = 1:numel(objA)
                    r(cnt) = tensorContainer(objA(e).ten + objB(e).ten);
                    cnt = cnt + 1;
                end
                r = reshape(r,size(objA));
            end
            
        end
        
        function [] = sum(obj,order)
            
        end
        
        % get the order of the tensor
        function [r] = order(obj)
            r = size(obj);
            r(r==1) = [];
            if isempty(r)
                r = 1;
            else
                 r = numel(r);
            end
        end
        
    end
    
    methods (Static)
        
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
        
        % apply operation along 
        function [varargout] = applyUniary(array,uniaryOp)
            % for each element
            for e = 1:numel(array)
                if nargout ~= 0
                    varargout{1}(e) = tensorContainer(uniaryOp(array(e).ten));
                else
                    uniaryOp(array(e).ten);
                end
            end
            % preserve the order
            if nargout ~= 0
                varargout{1} = reshape(varargout{1},size(array));
            end
        end
        
        function [varargout] = applyBinary(binaryOperand1,binaryOp,binaryOperand2)
            if isa(binaryOperand1,'tensorContainer') && isa(binaryOperand2,'tensorContainer')
                if tensorContainer.orderCheck(binaryOperand1,binaryOperand2)
                    if tensorContainer.dimsCheck(binaryOperand1,binaryOperand2)
                        % for each element
                        for e = 1:numel(binaryOperand1)
                            varargout{1}(e) = tensorContainer(binaryOp(binaryOperand1(e).ten,binaryOperand2(e).ten));
                        end
                        % preserve the order
                        varargout{1} = reshape(varargout{1},size(binaryOperand1));
                    end
                end
            end
        end
        
        % perform accumulation operation along order
        function [result] = applyDotBinary(array,binaryOp,order)
            % pull the operation order to the front order
            [array,sz,newOrder] = tensorContainer.pullFront(array,order);
            % for each column - operate on vector
            for column = 1:size(array,2)
                % perform accumulation
                if size(array,1) > 1
                    accumulator = binaryOp(array(1,column),array(2,column));
                    for row = 2:(size(array,1)-1)
                        accumulator = binaryOp(accumulator,array(row+1,column));
                    end
                else
                    accumulator = array(1,column);
                end
                result(column) = tensorContainer(accumulator);
            end
            % reshape the tensor container
            result = reshape(result,[sz(2:end) 1]);
        end
        
        function [] = applyUniaryDotBinary()
            
        end
        
        function [] = applyBinaryDotBinary()
        
        end
        
        function [] = applyTensorBinary()
        end
    end
end





%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%
    % test the addition of the tensor container object
    %%%%%%%%%%%%%%%%%%%%%%%%
    a1 = tensorContainer(tensor(rand(3,2,3)));
    a2 = tensorContainer(tensor(rand(3,2,3)));
    a3 = tensorContainer(tensor(rand(3,2,3)));

    A = [a1;a2;a3];
    B = [a2;a2;a2];
    C = [A B];
    D = [B B];
    E = [B B B];

    r1 = tensorContainer.orderCheck(a1,a2);
    r2 = tensorContainer.orderCheck(a1,a3);
    r3 = tensorContainer.dimsCheck(a1,a3);
    r4 = tensorContainer.orderCheck(A,B);
    r5 = tensorContainer.orderCheck(C,D);
    r6 = tensorContainer.dimsCheck(C,D,[1 1]);
    r7 = tensorContainer.dimsCheck(C,D,[2 2]);
    r8 = tensorContainer.dimsCheck(C,D,[1 1]);
    r9 = tensorContainer.dimsCheck(C,D,[[1 1];[2 2]]);


    %%%%%%%%%%%%%%%%%%%%%%%%
    % apply uniary operator
    %%%%%%%%%%%%%%%%%%%%%%%%
    % test sum
    t1 = tensorContainer(tensor(rand(3,4,5)));
    t2 = tensorContainer(tensor(rand(1,2,3)));
    T = [t1;t2];
    R = T.applyUniary(str2func('sum'),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % test the apply uniary operation to each tensor in the tensor
    % test dotProduct as uniary op with fixed element
`   %%%%%%%%%%%%%%%%%%%%%%%%
    t1 = tensorContainer(tensor(rand(3,4,5)));
    t2 = tensorContainer(tensor(rand(3,2,3)));
    t3 = tensorContainer(tensor(rand(3,2,3)));
    T = [t1;t2;t3];
    op = tensor(rand(3,2),[1 1]);
    opFunc = @(X)tensor.dotProduct(X,op,[1 1]);
    R = T.applyUniary(opFunc);

    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % test uniary op as sum along 1
    %%%%%%%%%%%%%%%%%%%%%%%%
    opFunc = @(X)applyAlong(X,@(X,Y)sum(X,Y),1);
    T.applyUniary(opFunc);
    


    %%%%%%%%%%%%%%%%%%%%%%%%
    % test the addition of the tensor container object
    %%%%%%%%%%%%%%%%%%%%%%%%
    a1 = tensorContainer(tensor(rand(3,2,3)));
    a2 = tensorContainer(tensor(rand(3,2,3)));
    a3 = tensorContainer(tensor(rand(3,2,3)));
    A = [a1;a2;a3];
    b1 = tensorContainer(tensor(rand(3,2,3)));
    b2 = tensorContainer(tensor(rand(3,2,3)));
    b3 = tensorContainer(tensor(rand(3,2,3)));
    B = [b1;b2;b3];
    C = A + B;
    D = [A,A] + [B,C];
    G = cat(3,D,[A,C]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % try dot binary
    % use the sum operation- as cumsum
    %%%%%%%%%%%%%%%%%%%%%%%%
    binaryOp = @(X,Y)plus(X,Y);
    R = T.applyDotBinary(G,binaryOp,1);



    a1 = tensorContainer(tensor(rand(3,2,3)));
    a2 = tensorContainer(tensor(rand(3,2,3)));
    a3 = tensorContainer(tensor(rand(3,2,3)));
    A = [a1;a2;a3];

    meanOp = @(X)sum(X,1);
    U = tensorContainer.applyUniary(A,meanOp);
    subtractOp = @(X,Y)tensor.bsxfun(@minus,X,Y,[0 -1]);
    As = tensorContainer.applyBinary(A,subtractOp,U);














%}