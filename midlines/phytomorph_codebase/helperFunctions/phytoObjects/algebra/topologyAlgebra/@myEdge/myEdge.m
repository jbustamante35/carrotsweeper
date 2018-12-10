classdef myEdge < myHS_X & sizeConstraint & geometrical & sampleable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % my edge
    % one for one type - note that there is an added constraint of
    % max objects
    properties

    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myEdge(varargin)
            % super constructorss
            obj = obj@geometrical();
            obj = obj@myHS_X('myHS');
            obj = obj@sizeConstraint(2);
            %%%%%%%%%%%%%
            obj.view_props.props.LineStyle = '--';
            obj.view_props.type = 'segmentSet';
            obj.view_props.Color = 'r';
            %%%%%%%%%%%%%
            % assign the set - problem
            if nargin == 1
                setSource(obj,varargin{1}{1});
                setTarget(obj,varargin{1}{2});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set source
        function [] = setSource(obj,source)
            S(1).type = '{}';
            S(1).subs{1} = 1;
            obj = subsasgn(obj,S,source);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set source
        function [] = setTarget(obj,target)
            S(1).type = '{}';
            S(1).subs{1} = 2;
            obj = subsasgn(obj,S,target);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set source
        function [ret] = getSource(obj)
            S(1).type = '{}';
            S(1).subs{1} = 1;
            ret = subsref(obj,S);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set source
        function [ret] = getTarget(obj)
            S(1).type = '{}';
            S(1).subs{1} = 2;
            ret = subsref(obj,S);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % allow index
        function [b] = allowIndex(obj,S)
            b = S(1).subs{1} <= 2;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % add elements to allowable list     
        function [b] = isAllow(obj,ele)
            b = isAllow@typeConstraint(obj,ele);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload subsasgn
        function [ret] = subsasgn(obj,S,B)
            ret = obj;
            if obj.allowIndex(S)
                ret = subsasgn@myHS_X(obj,S,B);
            else
                fprintf(['object@index: ' num2str(S(1).subs{1}) ' was rejected\n']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view edge overload super class
        function [] = view(obj,h,frame,vProps)
            % set frame
            %if isempty(frame);frame = obj.bf;end
            data = pointSequence(obj,frame);
            data = [data(1:end-1,1);data(1:end-1,2)];
            lineSegment = myT(data);
            lineSegment.setView(obj.view_props);
            lineSegment.view(h,frame,vProps);
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample image at domain
        function [s] = sample(obj,image,domain)               
            % obtain representation of object as affine
            affine = rep(obj,[],'phytoAaffine');
            s = affine.sample(image,domain);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation
        function [r] = rep(obj,frame,type)
            % set frame
            if isempty(frame);frame = obj.bf;end
            % set type
            if isempty(type);type = obj.pref_representation;end
          
            % switch on rep type
            switch type
                case 'phytoApoint'
                    % point stack
                    data = pointSequence(obj,frame);
                    % normalize the point
                    r = phytoApoint(mean(data,2));       
                case 'phytoAaffine'
                    % point stack
                    data = pointSequence(obj,frame);
                    tVec = diff(data,1,2);
                    tVec(1:2) = flipud(tVec(1:2));
                    tVec = tVec / norm(tVec);
                    nVec = [tVec(2);-tVec(1);0];
                    l = length(obj);
                    tVec = .5*l*tVec;
                    %tVec = 100*tVec;
                    pt = mean(data,2);
                    initA = [nVec,tVec,pt];
                    r = phytoAaffine(initA);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation
        function [data] = pointSequence(obj,frame)
            % get the source and target nodes
            source = obj.getSource();
            target = obj.getTarget();
            % get first objects @ the nodes
            source = source{1};
            target = target{1};
            % get the reps for source and target
            pt_source = source.rep(frame,'phytoApoint');
            pt_target = target.rep(frame,'phytoApoint');
            % point sequence
            data = [pt_source.d',pt_target.d'];
            %data(1:2,1:2) = flipud(data(1:2,1:2));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % length
        function [l] = length(obj)
            data = pointSequence(obj,[]);
            dL = diff(data,1,2);
            l = sum(dL.*dL).^.5;
        end
    end
end

%{
S = myN();    
S{1} = myT(ones(2,1));

T = myN();    
T{1} = myT(ones(2,1));

E = myEdge();
E{3} = S;
E{1} = S;
E.getTarget()
E.getSource()
E.setTarget(T);


h = figure;hold on;
E.view(h);
%}




