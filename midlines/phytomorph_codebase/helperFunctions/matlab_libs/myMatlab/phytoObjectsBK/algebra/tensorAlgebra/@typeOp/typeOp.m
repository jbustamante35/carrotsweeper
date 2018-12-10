classdef typeOp < matlab.mixin.Copyable

    properties
                
        %%%% binary op pair
        bOP;

        %%%% vars for searching graph
        target;      % search type or search size -  property sought
        %target2;      % container type or else statement operand
        
        %%%% results from graph search
        r;      % last result
        lab;    % labels
        
        %%%% vars for keeping track of location during search
        n;      % current iteration level
        N;      % max level
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = typeOp(varargin)
            %%%%%%%%%%%%%%%%
            %%%% init search type
            obj.target{1} = '';
            obj.target{2} = '';
            %%%%%%%%%%%%%%%%
            %%%% init boolean switch to type
            % case1     : if object is element type then stack on return
            %           : if object is container type then open and investigate
            % case2     : if object is sized as constraint and if proper
            %               type then return
            %           : if object is not properly sized and is container,
            %               then open and investigate
            %%%%%%%%%%%%%%%%
            obj.bOP = cell(2,1);
            
            %obj.bOP{1} = @(A,B)typeOp.isaType(A,B);
            %obj.bOP{2} = @(A,B)typeOp.isaType(A,B);
            
            %obj.bOP{1} = @(A,B)typeOp.isSizeType(A,B);
            %obj.bOP{2} = @(A,B)typeOp.isaType(A,B);
            
            %%%%%%%%%%%%%%%%
            %%%% history
            obj.r = {};
            obj.lab = {};
            %%%%%%%%%%%%%%%%
            obj.N = inf;
            obj.n = 1;
            %{
            %%%%%%%%%%%%%%%%
            %%%% if init search type - init
            if nargin >= 1
                obj.target{1} = varargin{1};
            end
            %%%%%%%%%%%%%%%%
            %%%% if init container type - init
            if nargin == 2
                obj.target{2} = varargin{2};
            end
            %}
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set max recursive level
        function [] = setN(obj,N)
            obj.N = N;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % attach binary OP
        function [] = attachBO(obj,bOp,target,n)
            obj.bOP{n} = bOp;
            obj.target{n} = target;
        end
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % go
        function [] = op(obj,S)
            % iterate over each object in cell array
            for e = 1:S.numel()
                % current element
                ele = S{e};
                
                % if search type then store in result
                if obj.bOP{1}(ele,obj.target{1})

                    obj.r{end+1} = ele;
                    obj.lab{end+1} = ['e:' num2str(e) '--l:' num2str(obj.n)];

                    
                % if container type - recursive call
                elseif obj.bOP{2}(ele,obj.target{2})

                    % increment the level count on the tree
                    obj.n = obj.n + 1;
                    % if we are "allowed" to go deeper, then do so
                    if obj.n < obj.N
                        % call recursivly
                        obj.op(ele);
                    end
                    % decrement the current recurse level
                    obj.n = obj.n - 1;
                end

            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % clear results history
        function [] = clearH(obj)
            obj.r = {};
            obj.lab = {};
        end
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%
        % is a type
        function [ret] = isaType(A,B)
            ret = isa(A,B);
        end
        %%%%%%%%%%%%%%%%%
        % contains a type
        function [ret] = conaType(A,B)
            ret = A.type();
            ret = ~isempty(intersect(B,ret));
        end
        %%%%%%%%%%%%%%%%%
        % is size
        function [ret] = isSize(A,B)
            ret = isequal(numel(A), B);
        end
        %%%%%%%%%%%%%%%%%
        % not is size
        function [ret] = nisSize(A,B)
            ret = ~obj.isSize(A,B);
        end
        %%%%%%%%%%%%%%%%%
        % is single element - note hardwire to 1
        function [ret] = isSizeType(A,B)
            ret = 0;
            if typeOp.isSize(A,1)
                ret = typeOp.isaType(A{1},B);
            end
        end
    end
    
end


%{
    
%}