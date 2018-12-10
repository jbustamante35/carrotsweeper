classdef TA < handle
    
    
    properties
        gradeNames = {};
        tensorGrades = {};
        foldLog = {};
        toLog = false;
    end
    
    
    methods
        function [obj] = TA(varargin)
            if nargin ~= 0
                obj.tensorGrades = varargin{1};
                obj.gradeNames = varargin{2};
            end
        end
        
        function [] = toggleLog(obj)
            obj.toLog = ~obj.toLog;
        end
        
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'{}')
                try
                    if numel(S) == 1
                        for e = 1:nargout
                            varargout{e} = obj.getGrade(S(1).subs{1}(e));
                        end
                        %{
                        % meant for subindex such as T{1}
                        [varargout{1:nargout}] = obj.getGrade(S(1).subs{1});
                        %}
                    else
                        % meant for subindex such as T{1}(:,:,1)
                        [varargout{1:nargout}] = builtin('subsref',obj.tensorGrades{S(1).subs{1}},S(2:end));
                    end
                catch ME
                    getReport(ME)
                    fprintf(['please index tensor algrbra properly.\n']);
                end
            elseif strcmp(S(1).type,'.')
                try
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                    %{
                    if numel(S) == 1
                        % assume field text level indexing
                        gradeIndex = find(strcmp(obj.gradeNames,S(1).subs));
                        varargout{1} = obj.getGrade(gradeIndex);
                    else
                        builtin('subsref',obj,S)
                    end
                    %}
                catch
                    fprintf(['please index tensor algrbra properly.\n']);
                end
            elseif strcmp(S(1).type,'()') 
                  % meant for subindex such as T{1}
                  [varargout{1:nargout}] = obj.getGrade(S(1).subs{1});
            else
                fprintf(['please index tensor algrbra properly.\n'])
            end
            here = 1;
        end
        
        % return the number of grades
        function [g] = nGrade(obj)
            g = numel(obj.tensorGrades);
        end
        
        % return the number of trials
        function [n] = numelN(obj)
            n = ndims(obj.tensorGrades{1});
            n = size(obj.tensorGrades{1},n);
        end
        
        function [] = putGrade(obj,grade,name)
            obj.tensorGrades{end+1} = grade;
            if nargin == 3
                obj.gradeNames{numel(obj.tensorGrades)} = name;
            end
        end
        
        function [ret] = getGrade(obj,name)
            
            
            if ischar(name)
                name = find(strcmp(obj.gradeNames,name));
            end
            ret = obj.tensorGrades{name};
            
            %ret = TA(obj.tensorGrades(name),obj.gradeNames(name));
            %{
            if ~iscell(name)
                name = num2cell(name);
            end
            for e = 1:numel(name)
                if ischar(name{e})
                    name{e} = find(strcmp(obj.gradeNames,name{e}));
                end
                varargout{e} = TA(obj.tensorGrades(name{e}),obj.gradeNames(name{e}));
            end
            %}
        end
        
        % for trial index idx for fixed grade
        function [ret] = getGradeN(obj,grade,idx)
            order = ndims(obj.tensorGrades{grade});
            [idx] = obj.renderTrialIndex(order,idx);
            ret = obj.tensorGrades{grade}(idx{:});
        end
        
        % get the idxth trial
        function [ret] = getN(obj,idx)
            for grade = 1:numel(obj.tensorGrades)
                subT{grade} = obj.getGradeN(grade,idx);
                subN{grade} = obj.gradeNames{grade};
            end
            ret = TA(subT,subN);
        end
        
        function [obj] = subsasgn(obj,S,data)
            if numel(S) == 1
                try
                    if S(1).subs{1} > numel(obj.tensorGrades)
                        obj.putGrade(S(1).subs(1),data);
                    else
                        obj.tensorGrades{S(1).subs{1}} = data;
                    end
                catch
                     fprintf(['please assign tensor algrbra properly.\n']);
                end
            else
                obj.tensorGrades{S(1).subs{1}} = builtin('subsasgn',obj.tensorGrades{S(1).subs{1}},S(2:end),data);
            end
        end
        
        function [] = deleteGrade(obj,name)
            if ischar(name)
                name = find(strcmp(obj.gradeNames,name));
            end
            obj.tensorGrades(name) = [];
        end
        
        function [] = cat(grade,dim,data)
            obj{grade} = cat(dim,obj{grade},data);
        end
        
        function [yData] = applyFunc(obj,func,grade,storeName,idx)
            if nargin <= 4
                idx = 1:obj.numelN();
            end
            tData = obj.getGradeN(grade,1);
            tData = func(tData);
            sz = size(tData);
            yData = zeros(prod(sz),numel(idx));
            for e = 1:numel(idx)
                tData = obj.getGradeN(grade,idx(e));
                y = func(tData);
                yData(:,e) = y(:);
            end
            yData = reshape(yData,[sz numel(idx)]);
            if ~isempty(storeName)
                obj.putGrade(yData,storeName);
            end
        end
        
        function [] = createMappingNetworks()
        end
        
        function [] = reshape(obj,grade,newSZ)
            obj.tensorGrades{grade} = reshape(obj.tensorGrades{grade},newSZ);
        end
        
        function [] = permute(obj,grade,newO)
            obj.tensorGrades{grade} = permute(obj.tensorGrades{grade},newO);
        end
        
        
    end
    
    
    methods (Static)
        % apply dot without respect the the typing of the tensor
        function [r] = dot(a,b,d)
            n = ndims(a);
            na = 1:ndims(a);
            nb = 1:ndims(b);
            na = [d(1) setdiff(na,d(1))];
            nb = [d(2) setdiff(nb,d(2))];
            sza0 = size(a);
            szb0 = size(b);
            a = permute(a,na);
            b = permute(b,nb);
            sza1 = size(a);
            szb1 = size(b);
            a = reshape(a,[sza1(1) prod(sza1(2:end))]);
            b = reshape(b,[szb1(1) prod(szb1(2:end))]);
            r = mtimesx(a,'T',b);
            r = reshape(r,[sza1(2:end) szb1(2:end)]);
            %{
            p = 1:n;
            toS = p(d(1));
            p(d(1)) = p(end);
            p(end) = toS;
            r = permute(r,p);
            %}
        end
        
        function [r] = tensor(a,b)
            r = mtimesx(a(:),b(:),'T');
            r = reshape(r,[size(a) size(b)]);
        end
    end
    
    
    methods (Access=private)
        function [idx] = renderTrialIndex(obj,order,n)
            idx = repmat({':'},1,order);
            idx{end} = n;
        end
    end
end

%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % roll back logs
    clear T
    clear classes TA
    T = TA;
    T.putGrade(rand(1,2,3,4),'test')
    %T.permute(1,[2 1 3 4]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ret = T.getGradeN(1,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % test tensor product
    a = rand(2,3);
    b = rand(4,1);
    ab = TA.tensor(a,b);
    c = rand(4,2,6);
    ac = TA.tensor(a,c);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    clear T
    clear classes TA
    T = TA;
    TR = 20;
    T.putGrade(rand(10,10,TR),'test');
    T.putGrade(rand(2,10,50,TR),'test2');
    T.putGrade(rand(1,TR),'test3');
    T.putGrade({'a','b','c'},'test4');
    T.numelN();
    func = @(X)sum(X,1);
    r = T.applyFunc(func,1,'hello');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ret = T.getN(3:6);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    subT = T{1:2};
    subTV = T(1:2);
    func = @(T,G)TA.dot(T{1},G{2},[1 2]);
    R = func(T,T);
%}