classdef indexable < matlab.mixin.Copyable
    
    properties 
        indexDim;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % helpers - private
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get from tensor
        function [r] = get(obj,S)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get the data from the block
            % if S has index values            
            if ~isempty(S.subs)
                % generate index
                r = genIndexKey(obj,S);
                % get the object via the index
                r = loopUp(r);
                % operate to shape the objects
                r = postShape(r);
            % if there is no index than baseR == 0
            else
                r = obj.d;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % put to tensor
        function [r] = put(obj,S,A)
            if ~isempty(S.subs)
                % generate index
                r = genIndexKey(obj,S);
                % operate to shape the objects
                A = preShape(A);
                % insert 
                insert(A,r);
            % if there is no index than baseR == 0
            else
                r = obj.d;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate polynomial index
        function [s] = gen_pi(obj)
            s = size(indexDim)
            s = cumprod(s);
            s = [1 s];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate index
        function [v] = genIndexKey(obj,S)
            v = 0;
            if ~isempty(S.subs)
                %%% generate kronsecker product of sets of index
                % if : as index
                v = colon(obj,S(1).subs{end},numel(S(1).subs));
                % raster S to index
                for e = 1:numel(S(1).subs)-1
                    curV = colon(obj,S(1).subs{end-e},numel(S(1).subs)-e);
                    mag = ones(1,numel(curV));
                    v = kron(v,mag);
                    v = [repmat(curV,[1 size(v,2)/size(curV,2)]);v];
                end
                %%% generate kronsecker product of sets of index
                % generate poly
                s = gen_pi(obj);
                % get poly part - aka remove the last product of dims
                s = s(1:size(v,1));
                % generate ind
                v = s*(v-1)+1;
            end 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload end
        function [r] = end(obj,k)
            r = obj.h(end).x(k);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % overload colon
        function [v] = colon(obj,v,e)
            % if colon as index then stack 1:N
            % where N is the dims of the vector
            % and N is stored in the history of the 
            % reshape
            if isa(v,'char')
                if strcmp(v,':')
                    v = 1:obj.h(end).x(e);
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % helpers - private
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods (Abstract)
        [ret] = preShape(ret);
        [] = insert(ret);
        [ret] = lookUp(ind);
        [ret] = postShape(ret);
    end
end