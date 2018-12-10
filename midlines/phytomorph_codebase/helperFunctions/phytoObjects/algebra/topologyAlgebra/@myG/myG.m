classdef myG < viewable & sampleable
    % my node class
    properties
        E;  % edges
        N;  % nodes
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myG(varargin)
            % init nodes and edges
            obj.N = myN();  % node for nodes
            obj.E = myHS_X('myEdge');
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert node
        function [idx] = putNode(obj,n,idx)
            if nargin == 2
                idx = obj.N.numel() + 1;
            end
            obj.N{idx} = n;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert edge
        function [idx] = putEdge(obj,e,idx)
            if nargin == 2
                idx = obj.E.numel() + 1;
            end
            obj.E{idx} = e;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % bulk insert nodes
        function [obj] = bputNode(obj,s)
            for e = 1:numel(s)
                n = myN();
                n{1} = s{e};
                putNode(obj,n,e);
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert node
        function [ret] = getNode(obj,idx)
            ret = obj.N{idx};
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert edge
        function [ret] = getEdge(obj,idx)
            ret = obj.E{idx};
        end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view nodes
        function [] = viewNodes(obj,h)
            obj.N.view(h);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % view edges
        function [] = viewEdges(obj,h)
            obj.E.view(h);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % count edges
        function [e_sz] = edgeSize(obj)
            e_sz = numel(obj.E);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % count nodes
        function [n_sz] = nodeSize(obj)
            n_sz = numel(obj.N);
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample image at domain
        function [] = sampleViewer(obj,domain,h,vProps)
            if nargin == 3;vProps.nodes = [];vProps.edges = [];end
            %%%%%%%%%%%%%%%%%%
            % edge
            for e = 1:numel(obj.E)
                % obtain representation of edge object as affine
                affine = rep(obj.E{e},[],'phytoAaffine');
                affine = affine.d;
                affine(1:2,1:2) = flipud(affine(1:2,1:2));
                viewD = affine*domain.edge.view_rep;
                % view curve
                curve = phytoAcurve(viewD');
                curve.view(h,[],vProps.edges);
            end
            %%%%%%%%%%%%%%%%%%
            % node
            for n = 1:numesubsref(r,S)l(obj.N)
                % get the node
                nd = obj.N{n};
                % obtain representation of object as affine
                affine = rep(nd.getElement(1),[],'phytoAaffine');
                viewD = affine.d*domain.node.view_rep;
                % view curve
                curve = phytoAcurve(viewD');
                curve.view(h,[],vProps.edges);
            end
            
            
        end
        %}
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % sample image at domain
        function [s] = sample(obj,image,domain)
            % create return sample graph
            s = gsS();
            % sample the edges
            for e = 1:numel(obj.E)
                sE = obj.E{e}.sample(image,domain.edge);
                s.insertEdge(sE);
            end
            % sample the edges
            for e = 1:numel(obj.N)
                sN = obj.N{e};
                sN = sN.getElement(1);
                sN = sN.sample(image,domain.node);
                s.insertNode(sN);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % representation        
        function [r] = rep(obj,frame,type)
            % set frame
            if isempty(frame);frame = obj.bf;end            
            % switch on rep type
            switch type
                case 'phytoPoint'
                   
                case 'phytoAffine'
                   
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % vectorize
        function [v] = vectorize(obj,v)
           v = [v obj.N.vectorize(v)];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % call view
        function [] = view(obj,h,frame,vProps)            
            if nargin == 3;vProps.nodes = [];vProps.edges = [];end
            obj.N.view(h,frame,vProps.nodes);
            obj.E.view(h,frame,vProps.edges);
        end
        function [r] = subsref(obj,S)
            % if . type
            if strcmp(S(1).type,'.')
                if any(strcmp(methods(obj),S(1).subs))
                    try r = builtin('subsref',obj,S(1:2));catch ME;builtin('subsref',obj,S(1:2));end
                    S(1:2) = [];
                else
                    try r = builtin('subsref',obj,S(1));catch ME;builtin('subsref',obj,S(1));end
                    S(1) = [];
                end
            % if () type    
            elseif strcmp(S(1).type,'()')
                try r = builtin('subsref',obj.S,S(1));catch ME;builtin('subsref',obj.S,S(1));end
                S(1) = [];
            % if {} type
            elseif strcmp(S(1).type,'{}')
                try r = builtin('subsref',obj.S,S(1));catch ME;builtin('subsref',obj.S,S(1));end
                S(1) = [];
            end
            % recursive call(s)
            if numel(S) >= 1
               % call to the next index level
               try r = subsref(r,S);
               catch ME;
                    try subsref(r,S);
                    catch ME% call to the next index level
                        try r = builtin('subsref',r,S);
                        catch ME;
                            builtin('subsref',r,S);
                        end
                    end
               end
            end
        end
    end
        
    methods (Static)
        %%%%%%%%%
        %%% Function parse_inputs
        %%%%%%%%%
        function [extraArgs, msg] = parse_inputs(v)
            extraArgs = [];
            msg = '';
            try
                if mod(numel(v),2) == 1
                    msg = 'Input ill formed.\n';
                else
                    % loop over other arguments
                    for e = 1:(numel(v)/2)
                        prop = (e-1)*2 + 1;
                        value = prop + 1;
                        extraArgs.(v{prop}) = v{value};
                    end
                end
            catch
                
            end
        end
    end
    
        
end

%{
%%%%%%%%%%%%%%%%
% graph
G = myG();

%%%%%%%%%%%%%%%%
% generate X rand nodes
dV = phytoPoint([1 1]');
N = 100;
for e = 1:N
    tN = myN();
    vl = rand(2,1);
    pt = phytoPoint(vl);
    %{
    E = eye(3);
    E(1:2,3) = vl;
    pt = phytoAffine(E);
    %}
    %pt.displace(dV);
    tN{1} = pt;
    G.putNode(tN);
end

%%%%%%%%%%%%%%%%
% generate M connections
M = 10;
for e = 1:M
    r1 = round(1 + (N-1).*rand(1));
    r2 = round(1 + (N-1).*rand(1));
    tE = myEdge();
    tE.setSource(G.getNode(r1));
    tE.setTarget(G.getNode(r2));
    G.putEdge(tE);
    SL(e) = r1;
end


h = figure;hold on;
G.view(h,[]);


h = figure;
p = @(t).2*[cos(t);sin(t)];
TH = linspace(-pi,pi,100);
for e = 1:numel(TH)
    dV = phytoPoint(p(TH(e)));
    dV.normalize(1);
    n = G.getNode(SL(1));
    pt = n{1};
    pt.displace(dV);
    hold on
    G.view(h,[]);
    drawnow
end


%}




