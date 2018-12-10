classdef tE < matlab.mixin.Heterogeneous
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this base object is a logical vector
    % the two states are:
    % false = co-vector
    % true = vector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        % type vector
        byte = false;
    end

    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = tE(initType,length)
            if nargin == 0
                initType= 0;
            end
            if issparse(initType)
                obj.byte = initType;
            else
                if nargin < 2
                    length = 1;
                end
                obj.byte = logical(sparse(length,numel(initType)));
                for e = 1:size(initType,2)
                    obj.byte(:,e) = logical(initType(e));
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % two states are equal iff
        % they are element-wise unequal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [r] = eq(objA,objB)
            try
                % for each object in the array
                for a = 1:numel(objA)
                    % test if all objects in the byte array 
                    % are equal
                    r(a) = all(all(objA(a).byte == ~objB(a).byte,2),1);
                end
                % all objects in array are equal
                r = full(all(r));
            catch ME
                fprintf(['error generated during type check at [tE] level.\n']);
            end
        end
        
        
        function [str] = toString(obj)
            poly = 2.^(0:size(obj.byte,2)-1);
            type = full(obj.byte)*poly';
            UQ = unique(type);
            str = ['(' num2str(size(obj.byte,2)) ':('];
            while ~isempty(type)
                msk = type(1) == type;
                msk = imfill(msk,1);
                v = type(1);
                l = sum(msk);
                str = [str num2str(v) '-' num2str(l) ')('];
                type(find(msk)) = [];
            end
            str(end) = [];
            str = [str ':' num2str(size(obj.byte,1))];
            str = [str ')'];
        end
        
        
        %{
        function [varargout] = glue(objA,objB)
            newByte = [];
            for a = 1:numel(objA)
                newByte = [newByte objA(a).byte];
            end
            if nargin == 2
                for b = 1:numel(objB)
                    newByte = [newByte objB(b).byte];
                end
            end
            varargout{1} = tE(newByte);
        end
        %}
        
        function [r] = fold(objA,objB)
            newA = repmat(objA.byte,[size(objB.byte,1) 1]);
            newB = repmat(objB.byte,[size(objA.byte,1) 1]);
            r = tE([newA newB]);
        end
        
        %
        function [r] = directSum(objA,objB)
            r = tE([objA.byte;objB.byte]);
        end
        
        %{
        function [] = swap(obj,p)
            for e = 1:numel(obj)
                obj(e).byte = obj(e).byte(p);
            end
        end
        %}
        
   end
end

%{
   
    

%}