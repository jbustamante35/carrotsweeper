classdef TE < matlab.mixin.Heterogeneous
   properties
       Byte;
   end
   
   methods
        function [obj] = TE(iByte)
            if nargin == 0
                iByte = tE(0);
            end
            obj.Byte = iByte;
        end
        
         
        function [r] = toString(obj)
            %r = '<';
            r = '';
            for e = 1:numel(obj)
                r = [r strrep(toString(obj(e).Byte),'  ','') '|'];
            end
            r(end) = [];
            %r = [r '|>'];
        end
        
        function [r] = eq(objA,objB)
            r = all(objA.Byte == ~objB.Byte);
        end
        
        
        function [r] = ByteSize(objA)
            r = zeros(numel(objA),1);
            for e = 1:numel(objA)
                r(e) = size(objA(1).Byte,2);
            end
        end
        
        % ask for the order from the signature
        function [r] = order(obj)
            r = size(obj);
            r(r==1) = [];
            r = numel(r);
        end
        
        % remove an order from the signature
        function [obj] = removeOrder(obj,order)
            obj(order) = [];
            %{
            % removed for verbose
            for e = 1:numel(obj)
                obj(e).Byte(order) = [];
            end
            r = repmat({':'},[1 obj.order]);
            r{order} = 1;
            obj = squeeze(obj(r{:}));
            %}
        end
        
        function [r] = getOrderType(obj,order)
            for e = 1:numel(obj)
                r(e) = obj(e).Byte(order);
            end
            r = reshape(r,size(obj));
        end
        
        % 
        function [varargout] = glue(objA,objB)
            szA = ByteSize(objA);
            szB = ByteSize(objB);
            cnt = 1;
            %newByte(sum(szA) + sum(szB)) = tE();
            newByte = [objA.Byte objB.Byte];
            %{
            for a = 1:numel(objA)
                    try
                    newByte(cnt) = objA(a).Byte;
                    cnt = cnt + 1;
                    catch ME
                        ME
                    end
                %newByte = [newByte objA(a).Byte];
            end
            for b = 1:numel(objB)
                newByte(cnt) = objA(a).Byte;
                cnt = cnt + 1;
                %newByte = [newByte objB(b).Byte];
            end
            %}
            varargout{1} = TE(newByte);
        end
        
        % swap the order of the signature
        function [obj] = swap(obj,p)
            obj = obj(p);
            %{
            for e = 1:numel(obj)
                obj(e).Byte = obj(e).Byte(p);
            end
            %}
        end
        
        % permute the signature for the element in the signature
        % and the order of the signature tensor
        function [] = permute(obj,order)
            %{
            % used for verbose signature
            for e = 1:numel(obj)
                obj(e).Byte = obj(e).Byte(order);
            end
            %}
            builtin('permute',obj,order);
        end
        % fold = permute and reshape by fixed values
        function [obj] = fold(obj,foldDirections)
            
            for f = 1:numel(foldDirections)
                %[objA.Byte objB.Byte]
                tmp = [obj(foldDirections{f}).Byte];
                newB(f) = TE(tmp.glue([]));
            end
            obj = newB;
            
            %{
            pVec = [];
            %oldSZ = (obj);
            newSZ = [];
            newF = {};
            % the "old" order vector is always the "natural"
            oldF = 1:numel(obj);
            % 
            for e = 1:numel(foldDirections)
                % create the permute vector
                pVec = [pVec foldDirections{e}];
                % create the new size vector
                newSZ = [newSZ prod(oldSZ(foldDirections{e}))];
                % create the new fold
                newF{e} = [oldF(foldDirections{e})];
            end
            % check for newSZ to have two elements
            if numel(newSZ) == 1
                newSZ = [newSZ 1];
            end
            % check to see if any dims are missed
            dF = setdiff(min(pVec):max(pVec),1:numel(oldSZ));
            % if the permutation is valid
            if isempty(dF)
               
                % call glue for the new signature
                
                for f = 1:numel(newF)
                    %[objA.Byte objB.Byte]
                    tmp = [obj(newF{f}).Byte];
                    newB(f) = tmp.glue([]);
                end
                obj = TE(newB);
                % 
                %obj = obj.swap(pVec);
                % reshape the signature 
                %obj = reshape(obj,newSZ);
            else
                fprintf(['The fold operation was not valid.\n'])
            end
            %}
        end
        
        %{
        % used for verbose signature - removed
        function [varargout] = mtimes(objA,objB)
            varargout{1}(numel(objB)*numel(objA)) = TE();
            cnt = 1;
            for b = 1:numel(objB)
                for a = 1:numel(objA)
                    f = [objA(a).Byte objB(b).Byte];
                    %varargout{1}(cnt) = glue(objA(a),objB(b));
                    %varargout{1}(cnt) = TE([objA(a).Byte objB(b).Byte]);
                    cnt = cnt + 1;
                end
                b
            end
            varargout{1} = squeeze(reshape(varargout{1},[size(objA) size(objB)]));
        end
        %}
        
        
        function [K] = tensorProduct(objA,objB)
            
        end
        
   end
end

%{


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % unit testing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create 
    t1 = tE([1]);
    t2 = tE([1]);
    t3 = tE([1]);
    
    t4 = tE([0]);
    t5 = tE([0]);



    T1 = TE([t1 t4]);
    T2 = TE([t2 t4]);
    T3 = TE([t3 t4]);

    T4 = TE([t1 t5]);
    T5 = TE([t2 t5]);
    T6 = TE([t3 t5]);

    T = [[T1;T2;T3],[T4;T5;T6]];

    %T.permute([2 1]);
    %T = T.fold({[1 2]});


    % test multiply
    R = T4*T6

    q1 = [T1;T2;T3];


    q2 = q1*T1;


    q1 = TE([1]);
    q2 = TE([1]);
    qT = [q1 q2];










%}