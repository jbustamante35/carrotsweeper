classdef phytoAbasisVector < myTB
    % class for set of basis vectors
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAbasisVector(varargin)
            obj@myTB();
            if (nargin == 1)
                % set the base rank, fibre rank, shapeFlag
                setAllRank(obj,0,ndims(varargin{1}));
                setData(obj,varargin{3});
                obj.sF = 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of vectors
        function [num] = vectorCount(obj)
            num = size(obj.d,2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dimension of vectors
        function [num] = vectorDim(obj)
            num = size(obj.d,1);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % tensorProduct
        function [] = tensor(obj,A)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn - assign bassis vectors
        function [obj] = subsasgn(obj,S,B)
            if strcmp(S(1).type,'.')                
                subsasgn@myTB(obj,S,B);
            else
                try
                    obj.d(:,S(1).subs{1}) = B;
                catch ME
                    ME
                end
            end
        end
        
    end
end