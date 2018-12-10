classdef tensor < handle
    properties
        M;
        signature = TE(tE());
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
                        %obj.signature = TE([]);
                        for idx = 1:typedIndexNumber
                            %vec(sM{idx}) = TE();
                            %obj.signature = obj.signature*vec;
                            obj.signature(idx) = TE(tE());
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
                        if numel(sig) == 1
                            sig = repmat(sig,[1 typedIndexNumber]);
                        end
                        %{
                        obj.signature = TE([]);
                        
                        for idx = 1:typedIndexNumber
                            vec(1:sM{idx}) = TE(tE(sig(idx)));
                            obj.signature = obj.signature*vec;
                        end
                        %}
                        for idx = 1:typedIndexNumber
                            %vec(sM{idx}) = TE();
                            %obj.signature = obj.signature*vec;
                            obj.signature(idx) = TE(tE(sig(idx)));
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
        
        function [r] = sig(obj,idx)
            if nargin == 1
                idx = 1:numel(obj.signature);
            end
            r = [];
            strg = '';
            for e = idx
                r = [r obj.signature(e).Byte];
                %{
                for t = 1:numel(tmp)
                    r{e} = tmp(t).byte;
                end
                %}
            end
        end
        
        function [r] = sig2str(obj)
            r = toString(obj.signature);
        end
        
        function [] = permute(obj,order)
            opString = ['<-w(' strrep(num2str(order),'  ',',') ')|'];
            preStr = obj.sig2str;
            % apply the permutation for the tensor
            obj.M = permute(obj.M,order);
            % apply the permatation for the signature
            obj.signature = obj.signature.swap(order);
            postStr = obj.sig2str;
            fprintf([postStr '=' opString preStr '|>.\n']);
        end
        
        % need to fix -might not need to call
        % leave empty  - only fill in if the
        % need for reshape is needed and folding does not work
        function [] = reshape(obj,newSZ)
           obj.M = reshape(obj.M,newSZ);
        end
        
        function [] = fold(obj,foldDirections)
            r = '{';
            for l = 1:numel(foldDirections)
                r = [r '[' strrep(num2str(foldDirections{l}),'  ',',') '],'];
            end
            r(end) = [];
            r = [r '}'];
            opString = ['<-f(' r ')|'];
            preStr = obj.sig2str();
            
            
            
            %obj.signature = obj.signature.fold(size(obj),foldDirections);
            naturalOrder = 1:obj.order();
            pVec = [];
            oldSZ = obj.size();
            newSZ = [];
            for e = 1:numel(foldDirections)
                pVec = [pVec foldDirections{e}];
                newSZ = [newSZ prod(oldSZ(foldDirections{e}))];
                %sigDirections{e} = 
            end
            newOrder = pVec(pVec);
            str = 1;
            for e = 1:numel(foldDirections)
                stp = (str+numel(foldDirections{e})-1);
                sigDirections{e} = str:stp;
                str = stp + 1;
            end
            dF = setdiff(min(pVec):max(pVec),1:numel(oldSZ));
            if isempty(dF)
                obj.permute(pVec);
                obj.reshape(newSZ);
            end
            obj.signature = obj.signature.fold(sigDirections);
            
            postStr = obj.sig2str();
            fprintf([postStr '=' opString preStr '|>.\n']);
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
        % should re-write for general pull to "front"
        % prolly should write push back again
        function [unFold,focusedFolded,vectorized,newSignature] = pullSignature(a,d)


        end
        
        
        
         
        
    end
    
    methods (Static)
        
        function [c] = typeCheck(a,b,d)
            c = a.sig(d(1)) == b.sig(d(2));
            
            %{
            % removed for verbose
            newSigA = [];
            newSigB = [];
            Aorder = order(a);
            Border = order(b);
            subA = a.signature.getOrderType(d(1));
            subB = b.signature.getOrderType(d(2));
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
            %}
        end
        
        function [c] = dimsCheck(a,b,d)
            szA = size(a.M);
            szB = size(b.M);
            c = szA(d(1)) == szB(d(2));
        end
        
        function [r] = dotProduct(a,b,d)
            opString = ['||(' num2str(d(1)) ')x(' num2str(d(2)) ')||'];
            if isa(a,'tensor') && isa(b,'tensor')
                if tensor.dimsCheck(a,b,d)
                    [typeC] = tensor.typeCheck(a,b,d);
                    if typeC
                        astr = a.sig2str();
                        bstr = b.sig2str();
                        at = a.signature;
                        bt = b.signature;
                        a = a.M;
                        b = b.M;
                        r = TA.dot(a,b,d);
                        at = at.removeOrder(d(1));
                        bt = bt.removeOrder(d(2));
                        % removed for verbose signature
                        %newt = at*bt;
                        newt = [at bt];
                        r = tensor(r,newt);
                    end

                end
                fprintf([r.sig2str() '= <|' astr opString bstr '|>.\n']);
            end

            
        end
        
        function [r] = tensorProduct(a,b)
            opString = ['||*||'];
            astr = a.sig2str();
            bstr = b.sig2str();
            r = TA.tensor(a.M,b.M);
            % removed for verbose
            %rSignature = a.signature*b.signature;
            rSignature = [a.signature b.signature];
            r = tensor(r,rSignature);
            fprintf([r.sig2str() '= <|' astr opString bstr '|>.\n']);
        end
    end
end
%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%
    % constructor tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % create a basic tensor with co-low signature
    t1 = tensor(rand(3,4,5));
    % create a basic tensor with high signature
    t2 = tensor(rand(3,1,2),1);
    % create a basic tensor with high signature
    t3 = tensor(rand(3,1),1);
    t4 = tensor(rand(1,3),1);
    % create a tensor with defined signature
    t5 = tensor(rand(3,4,5),[0 1 1]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % dot product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    % test the dot product along the dims
    r = tensor.dotProduct(t1,t3,[1 1]);
   
    r = tensor.dotProduct(t1,t2,[1 1]);
    r = tensor.dotProduct(t1,t5,[3 3]);

    %%%%%%%%%%%%%%%%%%%%%%%%
    % tensor product tests
    %%%%%%%%%%%%%%%%%%%%%%%%
    s = tensor.tensorProduct(t1,t3);
    q = tensor.tensorProduct(t1,t2);
    

    %%%%%%%%%%%%%%%%%%%%%%%%
    % permute tests  note that the signature needs to change too
    %%%%%%%%%%%%%%%%%%%%%%%%
    q.permute([1 2 4 3 5]);


    fT = tensor(rand(3,4,5,6,7,8),[1 1 1 0 0 1]);
    fT.fold({3 [4 5] [2 1 6]});
    
    t5 = multiLinearMap(rand(3,4,5),[0 1 1]);
    u2 = t5.getExpectedVector(2);
    t5.permute([2 1 3]);
    r1 = t5(u2)
    r2 = squeeze(mean(t5.M,1));





    q1 = tensor(rand(3,1),0);
    q2 = tensor(rand(3,1),1);
    OP = tensor(rand(3,3),[1 0]);
    r = <q1|OP|q2>;

%}






