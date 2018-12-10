classdef phytoAbasis < myTB
    methods
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAbasis(varargin)
            obj@myTB();
            if (nargin == 1)
                % set the base rank, fibre rank, shapeFlag
                setAllRank(obj,0,ndims(varargin{1}));
                setData(obj,varargin{3});
                obj.sF = 1;
            end
            if (nargin == 3)
                
            end
        end
        
        
        % overload subsasgn
        function [obj] = subsasgn(obj,S,B)
            if strcmp(S(1).type,'.')                
                subsasgn@myT(obj,S,B);                
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