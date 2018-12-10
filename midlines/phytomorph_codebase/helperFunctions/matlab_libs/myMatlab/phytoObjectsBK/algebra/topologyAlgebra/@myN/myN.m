classdef myN < myHS_X
    % my node class - each node can bundle many objects
    % default type is tensor - myT
    properties
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myN(varargin)
            % super constructor: a set with viewable objects
            % able to contain a container of self-type
            obj = obj@myHS_X({'viewable','myN'});
            %%%%%%%%%%%%%
            % init defaults
            obj.S = {};            
            %%%%%%%%%%%%%
            % init default views
            obj.view_props.props.MarkerType = '.';
            obj.view_props.props.MarkerColor = 'r';
            obj.view_props.type = 'phytoPoint';
            %%%%%%%%%%%%%
            % init default search
            obj.extractFunction = typeOp();
            % init op1 and attach
            op1 = @(A,B)typeOp.isSizeType(A,B);
            target1 = 'myT';
            obj.extractFunction.attachBO(op1,target1,1);
            % init op2 and attach
            op2 = @(A,B)typeOp.isaType(A,B);
            target2 = 'myN';
            obj.extractFunction.attachBO(op2,target2,2); 
            %%%%%%%%%%%%%
            % assign the set - problem
            if nargin >= 1
                obj.S = varargin{1};
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % default view implementation
        function [] = view(obj,h,frame,vProps)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % persist any props through function call            
            %userDef = viewable.parseView(varargin,obj.view_props.props);
            if nargin == 3;vProps = [];end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % harvest first level nodes
            obj.extractFunction.setN(1); 
            T = obj.harvest(obj.extractFunction);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % view each node - each node (default myT, must have view defined)
            for e = 1:numel(T)
                t = T{e};
                t = t.getElement(1);
                %t.setView(obj.view_props);
                t.view(h,frame,vProps);
            end
        end
    end
end

%{

S = myN();    
S{1} = 'h';
S{1} = myT(zeros(2,1));
S{2} = myT(ones(2,1));
%h = figure;S.view(h);


%}




