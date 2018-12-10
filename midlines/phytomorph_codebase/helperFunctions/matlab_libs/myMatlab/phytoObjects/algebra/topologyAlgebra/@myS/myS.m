classdef myS
    % set of objects - heterType
    % if set of sets are allowed,
    properties
        S = {};         % set        
        repFunc;        % default rep function        
        props = [];     % default view props
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = myS(varargin)
            obj.S = {};
            obj.repFunc = [];            
            if nargin == 1
                obj.S = varargin{1};
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % harvest
        function [r] = harvest(obj,op)
            % eval harvest over graph
            op.op(obj);
            % return results
            r = op.r;
            % clear res
            op.clearH();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % numel
        function [] = repFunction(obj,fun)            
            obj.repFunc = fun;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % numel
        function [r] = s(obj,varargin)
             r = prod(numel(obj.S));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload vertcat
        function [ret] = vertcat(varargin)
            ret = {};
            for e = 1:numel(varargin)                
                ret{e} = varargin{e};
            end
            ret = myS(ret);
        end      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload horzcat
        function [ret] =  horzcat(varargin)
            ret = {};
            for e = 1:numel(varargin)                
                ret{e} = varargin{e};
            end
            ret = myS(ret);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsref
        function [r] = subsref(obj,S)
            if strcmp(S(1).type,'.')
                if nargout == 0
                    builtin('subsref',obj,S);
                else
                    r = builtin('subsref',obj,S);
                end
            elseif strcmp(S(1).type,'()')
                r = builtin('subsref',obj.S,S); 
            elseif strcmp(S(1).type,'{}')
                e = builtin('subsref',obj.S,S);
                r = e;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [ret] = subsasgn(obj,S,B)
            if strcmp(S(1).type,'.')
                ret = builtin('subsasgn',obj,S,B);
            elseif strcmp(S(1).type,'()')
                obj.S = builtin('subsasgn',obj.S,S,B);
                ret = obj;
            elseif strcmp(S(1).type,'{}')
                obj.S = builtin('subsasgn',obj.S,S,B);                
                ret = obj;
            end
        end
        
    end
    
    methods (Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obtain integrated representation
        function [T] = rep(obj)
            T = obj.harvest(obj.repFunc);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view on figure h
        function [] = view(obj,h)
            % harvest first level nodes
            obj.repFunc.setN(1); 
            T = obj.harvest(obj.repFunc);
      
        end
    end
    
    methods (Static)
        function [props] = parseView(v)
            for e = 1:numel(v)/2
                strt = (e-1)*2+1;
                ed = str + 1;
                props.(v{strt}) = v{ed};
            end
        end
    end

end
%{
    % example with myS
    S = myS();    
    S{1} = 'h';S{2} = 1;S{3} = ['e',S{1}];S{4} = ['e',S{3}];
    cf = typeOp('char','myS');
    S.harvest(cf);
    cf.clearH();
    cf.setN(1);
    S.harvest(cf);

 
    S{1} = 'h';S{2} = 1;S{3} = ({'e',S{1}});S{4} = ({S{3},'w'});
    cO = typeOp('char','cell');
    cO.op(S);

%}

