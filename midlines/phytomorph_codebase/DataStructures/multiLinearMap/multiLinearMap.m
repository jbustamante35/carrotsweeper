classdef multiLinearMap < tensor
    properties
       
    end
    
    methods 
        
        function [obj] = multiLinearMap(varargin)
           obj = obj@tensor(varargin{:});
        end
        
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'{}')
                try 
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                catch ME
                    %getReport(ME)
                    fprintf(['error in index-{} for multilinear map.\n']);
                end
            elseif strcmp(S(1).type,'.')
                try
                    [varargout{1:nargout}] = builtin('subsref',obj,S);
                    
                catch
                    fprintf(['error in index-. for multilinear map.\n']);
                end
            elseif strcmp(S(1).type,'()')
                try
                    
                    %{
                    newM = obj.M;
                    % loop over the inputs
                    for e = 1:numel(S(1).subs)
                        if ~isempty(S(1).subs{e})
                            %{
                            sz = size(S(1).subs{e});
                            [v,d] = max(sz);
                            %}
                            sz = size(S(1).subs{e});
                            d = 1;
                            if sz(1) == size(obj.M,e)
                                %newM = TA.dot(newM,S(1).subs{e},[e d]);
                                newM = tensor.dot(newM,S(1).subs{e},[e d]);
                            else
                                error(['error: wrong vector size presented to the ' num2str(e)  'th input\n']);
                            end
                        else
                            
                        end
                    end
                    varargout{1} = multiLinearMap(newM);
                    %}
                    newM = obj;
                    % loop over the inputs
                    for e = 1:numel(S(1).subs)
                        if ~isempty(S(1).subs{e})
                            %{
                            sz = size(S(1).subs{e});
                            [v,d] = max(sz);
                            %}
                            sz = size(S(1).subs{e});
                            d = 1;
                            toOp = 1;
                            newM = tensor.dotProduct(newM,S(1).subs{e},[toOp d]);
                            %{
                            if sz(1) == size(obj.M,e)
                                %newM = TA.dot(newM,S(1).subs{e},[e d]);
                                newM = tensor.dot(newM,S(1).subs{e},[e d]);
                            else
                                error(['error: wrong vector size presented to the ' num2str(e)  'th input\n']);
                            end
                            %}
                        else
                            
                        end
                    end
                    varargout{1} = multiLinearMap(newM);
                catch ME
                    fprintf(['error in index-() for multilinear map.\n']);
                end
            else ME
                fprintf(['error in index-ANY for multilinear map.\n']);
            end
        end
        
    end
end


%{

    M = multiLinearMap(rand(3,3,3));
    v1 = tensor(rand(3,1),1);
    v2 = tensor(rand(3,1),1);
    v3 = rand(4,1);
    f = M(v1,v2);


    i1 = tensor([[0 0 1]' [1 0 0]'],1);
    %i1 = i1(:,1);
    v2 = tensor(rand(3,1),1);
    f = M(i1,v2);


    % test mean along dim1
    R = rand(4,20);
    U = mean(R,1);
    u = ones(1,size(R,1))/size(R,1);
    U2 = u*R;
    all(U == U2)
    RM = multiLinearMap(R);
    MV = tensor(u,1);
    MU = RM(MV);

    RM = multiLinearMap(R,[1 0]);
    RM.permute([2 1]);



    
    % example of index as function eval
    I = imread('/home/nate/Downloads/20170329n02_04.tif');
    I = multiLinearMap(I);

    


    
    

%}