classdef tensor < handle
    properties
        M;
        signature;
    end
    
    methods
          
        function [obj] = tensor(M,sig)
            if isa(M,'tensor')
                obj.M = M.M;
                obj.signature = M.signature;
            else


                M = squeeze(M);
                % transpose/permute - rows and column vectors are not things in this universe
                if size(M,1) == 1
                    M = permute(M,[2 1]);
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % calculate needed values
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nonTypedIndexNumber = sum(size(M) == 1);
                potentialTypedIndexNumber = ndims(M);
                typedIndexNumber = potentialTypedIndexNumber - nonTypedIndexNumber;
                sM = num2cell(size(M));
                if nargin == 1
                        % fun - empty set is my identity
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % assumed signature
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        obj.signature = TE([]);
                        for idx = 1:typedIndexNumber
                            vec(sM{idx}) = TE();
                            obj.signature = obj.signature*vec;
                        end
                        
                else
                    if ~isa(sig,'TE')
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % assumed signature
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % fun - empty set is my identity
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % assumed signature
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        obj.signature = TE([]);
                        if numel(sig) == 1
                            sig = repmat(sig,[1 typedIndexNumber]);
                        end
                        for idx = 1:typedIndexNumber
                            vec(1:sM{idx}) = TE(tE(sig(idx)));
                            obj.signature = obj.signature*vec;
                        end
                    else
                        obj.signature = sig; 
                    end

                end
                if ~isa(M,'double')
                    M = double(M);
                end
                obj.M = M;
            end
        end
        
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'{}')
                 fprintf(['error in index-{} for tensor.\n']);
            elseif strcmp(S(1).type,'.')
                try
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                catch ME
                    getReport(ME)
                    fprintf(['error in index-. for tensor.\n']);
                end
            elseif strcmp(S(1).type,'()')
                [varargout{1:nargout}] = builtin('subsref',obj.M,S);
            end
        end
        
        function [r] = size(obj)
            r = size(obj.M);
        end
        
        function [r] = order(obj)
            r = size(obj.M);
            r(r==1) = [];
            r = numel(r);
        end
        
        function [] = permute(obj,order)
            obj.M = permute(obj.M,order);
            obj.signature = obj.signature(:,order);
        end
        
        function [] = reshape(obj,newSZ)
            sig = full(obj.signature);
            oldSZ = size(obj.M);
            reshape(sig,[oldSZ size(sig,3)]);
        end
        
        function [] = fold(obj,foldDirections)
            pVec = [];
            oldSZ = obj.size();
            newSZ = [];
            for e = 1:numel(foldDirections)
                pVec = [pVec foldDirections{e}];
                newSZ = [newSZ prod(oldSZ(foldDirections{e}))];
            end
            dF = setdiff(min(pVec):max(pVec),1:numel(oldSZ));
            if isempty(dF)
                obj.permute(pVec);
                obj.reshape(newSZ);
            end
           
        end
        
        function [r] = sum(obj,order)
            r = tensor(sum(obj.M,order),obj.signature.removeOrder(order));
        end
        
        function [r] = prod(obj,order)
            r = tensor(prod(obj.M,order),obj.signature.removeOrder(order));
        end
        
        function [] = applyAlong(obj,func,toOpAlong)
            obj.M = func(obj.M,toOpAlong);
            obj.signature.removeOrder(order);
        end
        
        function [r] = or(objA,objB)
            test = 1;
        end
        
        function [r] = gt(objA,objB)
            test = 1;
        end
        
        function [r] = getExpectedVector(obj,d)
            sz = size(obj.M);
            v = ones(sz(d),1)/sz(d);
            [~,sig,~,~] = pullSignature(obj,d);
            r = tensor(v,sparse(~sig(:,1)));
        end
    end
    
    methods (Access = private)
         function [unFold,focusedFolded,vectorized,newSignature] = pullSignature(a,d)
            try

                 as = a.signature(:,d);
                 sz = size(a.M);
                 naturalIDX = 1:size(a.signature,2);
                 pullIDX = [d setdiff(naturalIDX,d)];
                 as = reshape(full(as),sz);

                 if numel(pullIDX) ~= 1
                    unFold = permute(as,pullIDX);
                 else
                    unFold = as;
                 end
                 newSZ = size(unFold);
                 focusedFolded = reshape(unFold,[newSZ(1) prod(newSZ(2:end))]);
                 vectorized = unFold(:);


                 if numel(pullIDX) > 1
                     for e = 2:numel(pullIDX)
                         nd = pullIDX(e);
                         as = a.signature(:,nd);
                         sz = size(a.M);
                         %naturalIDX = 1:size(a.signature,2);
                         %pullIDX = [nd setdiff(naturalIDX,nd)];
                         as = reshape(full(as),sz);


                         TMPunFold = permute(as,pullIDX);
                         TMPnewSZ = size(TMPunFold);
                         TMPfocusedFolded = reshape(TMPunFold,[TMPnewSZ(1) prod(TMPnewSZ(2:end))]);
                         TMPnewSignature = reshape(TMPfocusedFolded(1,:),[1 TMPnewSZ(2:end)]);
                         newSignature(:,e-1) = TMPnewSignature(:);
                     end
                 else
                     newSignature = sparse(1,1);
                 end
            catch ME
                getReport(ME)
            end


        end
        
    end
    
    methods (Static)
        
        function [c] = typeCheck(a,b,d)
            newSigA = [];
            newSigB = [];
            Aorder = order(a);
            Border = order(b);
            subA = a.signature.getOrderType(d(1));
            subB = b.signature.getOrderType(d(2));
            %{
            wholeA = [a.signature.Byte];
            wholeB = [b.signature.Byte];
            subA = reshape(wholeA(d(1):Aorder:end),size(a));
            subB = reshape(wholeB(d(1):Border:end),size(b));
            %}
            % pull dim to front for a,b
            newAorder = [d(1) setdiff(1:Aorder,d(1))];
            newBorder = [d(2) setdiff(1:Border,d(2))];
            if numel(newAorder) == 1
                newAorder = [newAorder 2];
            end
            if numel(newBorder) == 1
                newBorder = [newBorder 2];
            end
            subA = permute(subA,newAorder);
            subB = permute(subB,newBorder);
            newSZA = size(subA);
            newSZB = size(subB);
            subA = reshape(subA,[newSZA(1) prod(newSZA(2:end))]);
            subB = reshape(subB,[newSZB(1) prod(newSZB(2:end))]);
            % does the first element of each vector space match - pass for now
            c = subA(:,1) == subB(:,1);
            %r = bsxfun(@eq,subA,subB(:,1));
            %{
            [unFoldA,focusedFoldedA,vectorizedA,newSigA]  = pullSignature(a,d(1));
            [unFoldB,focusedFoldedB,vectorizedB,newSigB]  = pullSignature(b,d(2));
            %}
            %{
            % get the signature along the order d(x)
            as = a.signature(:,d(1));
            bs = b.signature(:,d(2));
            
            
            % get the natural size
            OszA = size(a.M);
            OszB = size(b.M);
           
            % pull the "operated-on" order/way to the front
            newDIM_order_a = [d(1) setdiff(1:numel(OszA),d(1))];
            newDIM_order_b = [d(2) setdiff(1:numel(OszB),d(2))];
            
            %%%%%%%%%%%%%%%%%%%%%%%%
            % make single-order/vector field of types
            %%%%%%%%%%%%%%%%%%%%%%%%
            % make full for type check
            as = reshape(full(as),OszA);
            % permute according to above
            as = permute(as,newDIM_order_a);
            % permute the size
            szA = size(as);
            
            
               
            %%%%%%%%%%%%%%%%%%%%%%%%
            % make single-order/vector field of types
            %%%%%%%%%%%%%%%%%%%%%%%%
            % make full for type check
            bs = reshape(full(bs),OszB);
            % permute according to above
            bs = permute(bs,newDIM_order_b);
            % permute the size
            szB = size(bs);
            
            % reshape for type check
            as = reshape(as,[szA(1) prod(szA(2:end))]);
            bs = reshape(bs,[szB(1) prod(szB(2:end))]);
            
            %}
            %{
            %%%%%%%%%%%%%%%%%%%%%%%%
            % homogenous requirement
            % the type must be the same orthogoal to the operation direction
            %%%%%%%%%%%%%%%%%%%%%%%%
            ca = all(bsxfun(@eq,focusedFoldedA,focusedFoldedA(:,1)),2);
            cb = all(bsxfun(@eq,focusedFoldedB,focusedFoldedB(:,1)),2);
            %}
            
            %{
            % set c to false
            c = false;
            % if the type matrix is homogenous
            if all(ca) & all(cb)
                % check the head of the list
                c = all(~focusedFoldedA(:,1) == focusedFoldedB(:,1));
                %{
                %{
                vecA = perms(newDIM_order_a(2:end));
                vecB = perms(newDIM_order_b(2:end));
                %}
               
                vecA = perms(2:numel(newDIM_order_a));
                vecB = perms(2:numel(newDIM_order_b));
                
                
                
                for e = 1:size(vecA,1)
                    newSigA = [newSigA reshape(as(1,:)',[1 prod(szA(vecA(e,:)))])'];
                end
               
                for e = 1:size(vecB,1)
                    newSigB = [newSigB reshape(bs(1,:)',[1 prod(szB(vecB(e,:)))])'];
                end
                %}
            end
            %c = ~strcmp(a.signature(d(1)),b.signature(d(2)));
            %}
        end
        
        function [c] = dimsCheck(a,b,d)
            szA = size(a.M);
            szB = size(b.M);
            c = szA(d(1)) == szB(d(2));
        end
        
        function [r] = dotProduct(a,b,d)
            if isa(a,'tensor') && isa(b,'tensor')
                if tensor.dimsCheck(a,b,d)
                    [typeC] = tensor.typeCheck(a,b,d);
                    if typeC

                        at = a.signature;
                        bt = b.signature;
                        a = a.M;
                        b = b.M;
                        r = TA.dot(a,b,d);
                        at = at.removeOrder(d(1));
                        bt = bt.removeOrder(d(2));
                        
                        %{
                        toBULKfromB = size(bt,1);
                        toBULKfromB = toBULKfromB*~isempty(bt);
                        if toBULKfromB == 0
                            toBULKfromB = 1;
                        end
                        
                        toBULKfromA = size(at,1);
                        toBULKfromA = toBULKfromA*~isempty(at);
                        if toBULKfromA == 0
                            toBULKfromA = 1;
                        end
                        %}
                        %{
                        newSigA_final = repmat(newSigA,[1*size(newSigB,1) 1]);
                        newSigB_final = repmat(newSigB,[1*size(newSigA,1) 1]);
                        if size(bt,2) == 0
                            newSigB_final = [];
                        end
                        newsig = [newSigA_final newSigB_final];
                        %}
                        newt = at*bt;
                        r = tensor(r,newt);
                    end

                end
            end

            
        end
        
        function [r] = tensorProduct(a,b)
            r = TA.tensor(a.M,b.M);
            rSignature = a.signature*b.signature;
            r = tensor(r,rSignature);
        end
    end
end
%{


    t1 = tensor(rand(3,4,5));
    t2 = tensor(rand(3,1,2),1);
    t3 = tensor(rand(3,1),1);
    t4 = tensor(rand(1,3),1);
    t5 = tensor(rand(3,4,5),[0 1 1]);

    %fT = tensor(rand(3,4,5,6,7,8));
    %fT.fold({3 [4 5] [2 1 6]});
    

    r = tensor.dotProduct(t1,t3,[1 1]);
    r = tensor.dotProduct(t1,t2,[1 1]);
    r = tensor.dotProduct(t1,t5,[3 3]);


    s = tensor.tensorProduct(t1,t3);
    q = tensor.tensorProduct(t1,t2);
    
    t5 = multiLinearMap(rand(3,4,5),[0 1 1]);
    u2 = t5.getExpectedVector(2);
    t5.permute([2 1 3]);
    r1 = t5(u2)
    r2 = squeeze(mean(t5.M,1));


    q.permute([1 2 4 3 5]);



    q1 = tensor(rand(3,1),0);
    q2 = tensor(rand(3,1),1);
    OP = tensor(rand(3,3),[1 0]);
    r = <q1|OP|q2>;

%}






