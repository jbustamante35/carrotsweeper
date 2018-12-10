classdef fo < matlab.mixin.Copyable%handle
    % function object
    properties
        ptDetector;
        key;
        moniker;
        func_handle;
        fmo_output_moniker;
        fmo_output_type;
    end
    
    methods
        function [obj] = fo(ptDetector,func_handle,moniker,fmo_output_moniker)
            
            
            obj.ptDetector = ptDetector; % make the point detector object available to the function object
            obj.func_handle = func_handle;
            obj.moniker = moniker;
            obj.key = obj.ptDetector.addFunction(obj);
            obj.fmo_output_moniker = fmo_output_moniker;
            
            
            
        end
        
        function [varargout] = run(obj,varargin)
            
            
            % generate the output key(s)
            if isa(varargin{1},'imageFile')
                [pth,nm,ext] = fileparts(varargin{1}.fileName);
            else
                nm = varargin{1}.getHeadNodeName();
            end
            
            
            % check if the output keys exist
            for e = 1:nargout
                outKeys{e} = ['hnn_' nm '_typeKey_' obj.fmo_output_type{e}];
                out(e) = exist([obj.ptDetector.oPath outKeys{e} '.mat']);
            end
            
            
            % if not all keys are computed then operator else load from key
            if ~all(out)
                % run the function
                varargout = obj.func_handle(varargin{:});
                % set the output type, persist, deposit
                for e = 1:numel(varargout)
                    varargout{e}.setTypeKey(obj.fmo_output_type{e});
                    varargout{e}.setptde(obj.ptDetector);
                    varargout{e}.setMoniker(obj.fmo_output_moniker{e});
                    % persist and deposit
                    varargout{e}.persist();
                    % store the feature map output
                    obj.ptDetector.depositFeatureMaps(varargout{e});
                    % store the ptDetector
                    obj.ptDetector.persist();
                end
            else
                for e = 1:numel(outKeys)
                    varargout{e} = obj.ptDetector.withdrawFeatureMap(outKeys{e});
                end
                fprintf(['function object computed previously \n']);
            end
        end
        
        function [varargout] = subsref(obj,s)
            flag = 1;
            if all(s(1).type == '()')
                if numel(s(1).subs) == 1
                    if isa(s(1).subs{1},'fmo')
                        flag = 0;
                        for e = 1:nargout
                            varargout{e} = fmo(obj.ptDetector,obj.fmo_output_moniker{e},'type');
                            obj.fmo_output_type{e} = varargout{e}.typeKey;
                        end
                    end
                end         
            end
            if flag
                if nargout == 0
                    builtin('subsref',obj,s);
                else
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                     %= {varargout};
                end
            end
        end
        
        function [] = persist(obj)
        
        end
        
        function [n] = getNargout(obj)
            n = numel(obj.fmo_output_type{e});
        end
    end
    
    methods (Static)
        function [] = renderClassType(conn,vf)
            import org.openrdf.model.vocabulary.RDF;
            import org.openrdf.model.vocabulary.RDFS;
            % make null
            cont = javaArray('org.openrdf.model.Resource',1);
            cont(1) = [];
            % make grand name space
            psNameSpace = 'http://example.org/productionSystem/';
            % make a feature map object
            fmoClass = vf.createURI(psNameSpace,'functionObject');
            conn.add(fmoClass,RDF.TYPE,RDFS.CLASS,cont);
        end
    end
end