classdef TE < matlab.mixin.Heterogeneous
   properties
       Byte;
   end
   
   methods
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = TE(iByte)
            if nargin == 0
                iByte = tE(0);
            end
            obj.Byte = iByte;
        end
        
         
        function [r] = toString(obj)
            r = '';
            for e = 1:numel(obj)
                r = [r strrep(toString(obj(e).Byte),'  ','') '|'];
            end
            r(end) = [];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the size of the tensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [r] = size(obj)
            r = size(obj.Byte);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ask for the order from the signature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [r] = order(obj)
            r = size(obj);
            r(r==1) = [];
            r = numel(r);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % remove an order from the signature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = removeOrder(obj,order)
            obj(order) = [];
        end
        
        function [r] = getOrderType(obj,order)
            for e = 1:numel(obj)
                r(e) = obj(e).Byte(order);
            end
            r = reshape(r,size(obj));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % swap the order of the signature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = swap(obj,p)
            obj = obj(p);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute the signature for the element in the signature
        % and the order of the signature tensor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [] = permute(obj,order)
            builtin('permute',obj,order);
        end
        
        % fold = permute and reshape by fixed values
        function [obj] = fold(obj,foldDirections)
            for f = 1:numel(foldDirections)
                tmp = [obj(foldDirections{f}).Byte];
                newB(f) = tmp.fold();
            end
            obj = newB;
        end
        
        
        function [K] = tensorProduct(objA,objB)
            
        end
        
   end
end

%{


    










%}