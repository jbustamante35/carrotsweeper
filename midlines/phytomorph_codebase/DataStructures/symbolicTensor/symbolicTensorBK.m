classdef symbolicTensor
    
    properties
        order;
        forwardMap;
        inverseMap;
    end
    
    methods
        function [obj] = symbolicTensor(map,infunc)
            order = size(map);
            
            uniMap = buildSymbolicTensor(order);
            uniMap = fillTensorComponents(uniMap,map,order);
            uniMap = renameTensorArgs(uniMap,order,['O' 'I']);
            
            if nargin == 2
                OS = sym(['O' '_%d'],[order(2) 1],'real');
                IS = sym(['I' '_%d'],[order(2) 1],'real');
                
                source = [OS;IS];
                target = [OS;infunc];
                uniMap = subs(uniMap,source,target);
            end
            
            obj.forwardMap = uniMap;
            obj.order = order;
            
            
        end
        
        function [varargout] = subsref(obj,S)
            if strcmp(S(1).type,'.')
                [varargout{1:nargout}] = builtin('subsref',obj,S);
            elseif strcmp(S(1).type,'()')
                
                
                if isa(S(1).subs{1},'smoothDomain')

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % make funcs on the fly
                    OS = eye(obj.order(1));
                    IS = sym(['I' '_%d'],[obj.order(2) 1],'real');
                    for s1 = 1:obj.order(1)
                        source = [];
                        for s2 = 1:obj.order(1)
                            source = [source;OS(s1,s2)];
                        end
                        source = [source;IS];

                        for e = 1:numel(source)
                            V{e} = source(e);
                        end
                        G{s1} = matlabFunction(obj.forwardMap(V{:}));
                    end

                    % make inputs
                    cnt = 1;
                    for ns = 1:numel(S(1).subs)
                        for s = 1:size(S(1).subs{ns}.M,3)
                            IN{cnt} = S(1).subs{ns}.M(:,:,s);
                            cnt = cnt + 1;
                        end
                    end

                    % eval
                    for e = 1:numel(G)
                        Y(:,:,e) = G{e}(IN{:});
                    end

                    varargout{1} = Y;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                elseif isa(S(1).subs{1},'symbolicTensor')
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    OS = sym(['O' '_%d'],[S(1).subs{1}.order(1) 1],'real');
                    IS = sym(['I' '_%d'],[obj.order(2) 1],'real');
                    TOT = [OS;IS];
                    map1 = S(1).subs{1}.forwardMap;
                    map2 = obj.forwardMap;
                    % make temp input and output
                    tmpOS = sym(['tmpO' '_%d'],[S(1).subs{1}.order(1) 1],'real');
                    tmpIS = sym(['tmpI' '_%d'],[obj.order(2) 1],'real');
                    tmpTOT1 = [tmpOS;IS];
                    tmpTOT2 = [OS;tmpIS];
                    map1(tmpTOT1) = subs(map1,TOT,tmpTOT1);
                    map2(tmpTOT2) = subs(map2,TOT,tmpTOT2);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    ostim = eye(S(1).subs{1}.order(1));
                    for o1 = 1:S(1).subs{1}.order(1)
                        stimVec = [];
                        for o2 = 1:S(1).subs{1}.order(1)
                            stimVec = [stimVec;ostim(1,o2)];
                        end
                        stimVec = [stimVec;IS];
                        IN = {};
                        for e = 1:numel(stimVec)
                            IN{e} = stimVec(e);
                        end
                        tmpMap(tmpOS(o1)) = map1(IN{:});
                        map2(TOT) = compose(map2,tmpMap,tmpIS(o1),tmpOS(o1));
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % make return
                    varargout{1} = symbolicTensor(map2);

                end
                
                
                
                
                
            end
        end
        
        
        
        
        
        
        
        
        
        
    end
end