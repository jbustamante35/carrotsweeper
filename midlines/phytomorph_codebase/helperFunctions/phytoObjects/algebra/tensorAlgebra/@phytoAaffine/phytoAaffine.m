classdef phytoAaffine < myTb %< phytoAgeo & sampleable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % my affine transformation
    % an affine transformation can be anywhere and anyway
    % it can be displaced and translated as an object
    % without changing its action
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = phytoAaffine(varargin)
            % super constructor - 0 rank-base and 2 rank-fibre
            obj = obj@myTb(0,2);
            % set default views
            obj.view_props.props.LineStyle = '--';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAaffine';
            % init point(s)
            if nargin == 1               
               obj.setData(varargin{1});
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view override
        function [] = view(obj,h,vProps)
            % view affine
            S.type = '()';
            S.subs = {};            
            uProps = obj.view_props;
            uProps.props = viewable.setProps(uProps.props,vProps);            
            tensorView(subsref(obj,S),uProps,h);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation        
        function [r] = rep(obj,frame,type)
            % set frame
            %if isempty(frame);frame = obj.bf;end
            % create eye
            d = eye(size(obj.d,1));
            % switch on rep type
            switch type
                case 'phytoApoint'
                    d = obj();
                    r = phytoApoint(d);                 
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample image at domain
        function [s] = sample(obj,image,domain)                                    
            tmp = domain.gen_parameters;
            I = zeros(tmp.sz);
            try
                % load para for q-reader
                para{1} = obj.d;
                para{2} = domain.d;
                sz = domain.gen_parameters;
                para{3} = sz.sz;
                I = myReader(image,'iatP',para);
            catch ME
                ME;
            end
            % return patch
            s = phytoApatch(I);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % distribute operator over the fibre and trial vectors
        function [] = distrib(obj,fibreI,trialI,op)
            % over each fibreI and each trialI
            for f = 1:numel(fibreI)
                for t = 1:numel(trialI)
                    tmp = getTrial(obj,trialI(t),fibreI(f));
                    tmp = op(tmp);
                    putTrial(obj,tmp,trialI(t),fibreI(f),0);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % angle
        function [data] = angle(obj)
           data = atan2(obj.d(1),obj.d(2));
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pull back
        function [id] = pullBack(d)
            id = inv(d);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % push forwards 
        function [ret] = pushForward(d,v)
            ret = d*v;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct from A point
        function [cur] = contructFromApoint(p)
            d = p();
            dA = eye(numel(d));
            dA(:,end) = d;
            cur = phytoAaffine(dA);
        end
    end
end

%{

    %%%%%%%%%%%%%%%%
    % graph
    G = myG();

    %%%%%%%%%%%%%%%%
    % generate X rand nodes
    N = 100;
    for e = 1:N
        tN = myN();
        pt = phytoPoint(rand(2,1));
        pt.normalize();
        G.putNode(tN);
    end

    ten = myT();



 %{
        
        %}
%}