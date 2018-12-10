classdef myTB < myHS_X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set of tensorBundles
    % first index will be baseTrial
    % 2-n index will be fibreTrial
    properties
        bR;             % base rank
        bI = 1;         % base index
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = myTB(varargin)
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
            addAllow(obj,'myTbS');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % need to set the base rank
            if nargin >= 1
                setBaseRank(obj,varargin{1});
            end
            % set the first data if one is provided
            % constructor can be: 
            %   1) a single myTbS
            %   2) a cell array of myTbS
            %   3) a single fibre
            if nargin >= 2
                
            end
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [ret] = subsasgn(obj,S,B)
            ret = [];
            try
                if checkBaseRank(obj,B)
                    ret = subsasgn@myHS_X(obj,S,B);
                else
                    fprintf(['object@class: ' class(B) ' was rejected\n']);
                end
            catch ME
                fprintf(['object@class: ' class(B) ' was rejected\n']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check for match
        function [b] = checkBaseRank(obj,fb)
            if ~isempty(obj.bR)
                b = (obj.bR == fb.getBaseRank());
            else
                fprintf(['need to set the base rank \n']);
            end
        end
    end

end