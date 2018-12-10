classdef myTbS < myHS_X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set of tensorBundles
    % first index will be baseTrial
    % 2-n index will be fibreTrial
    properties
        bR;             % base rank
        fI = 1;         % fibre index
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = myTbS(varargin)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % init calls as [baseRank,fibreRank]
            %   example - myTbS(1,1=3) is a curve in three space
            %                   a curve with the fibre structure of a 
            %                   three dim vector space
            %   example - myTbS(1,2 = [3 3]) is a curve with a second rank
            %                   tensor attach to the curve. the tensor is
            %                   R^3 tensor R^2
            %   eample - myTbS(2,1=3) is a surface in three space. the
            %                   along the surface is a vector fibre.
            % call super constructor
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj = obj@myHS_X();
            addAllow(obj,'myTb');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % need to set the base rank
            if nargin >= 1
                obj.bR = varargin{1};
                addBase(obj,varargin{2});
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % set the data for the base fibre
            if nargin >= 3
                obj.S{1}.setData(varargin{2});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add base fibre
        function [] = addBase(obj,fibreRank,baseRank)
            if nargin == 3;obj.bR = baseRank;end
            addFibre(obj,fibreRank,1);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add fibre
        function [] = addFibre(obj,fibreRank,index)
            obj.S{index} = myTb(obj.bR,fibreRank);
        end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set base rank
        function [r] = setBaseRank(obj,r)
            obj.bR = r;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get base rank
        function [r] = getBaseRank(obj)
            r = obj.bR;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % number of fibre spaces
        function [n] = nFibre(obj)
            n = size(obj.g0,obj.fI);
        end
        
    end

end