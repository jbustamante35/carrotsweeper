classdef featureFunctionBank < handle
    
    properties
        funcBank = [];
        monikerBank = [];
    end
    
    methods
        function [obj] = featureFunctionBank(varargin)
            if nargin == 1
                obj.funcBank = load(bankFile,'funcBank');
            end
        end
        
        function [key] = addFunction(obj,funcObject)
            funcstr = func2str(funcObject.func_handle);
            key = ['k_' num2str(string2hash(funcstr))];
            obj.funcBank.(key) = funcObject;
        end
        
        function [func] = getFunction(obj,key)
            func = obj.funcBank.(key);
        end
        
        function [] = removeFunction(obj,key)
            obj.funcBank.(key) = [];
        end
        
        function [] = persist(obj,location)
            funcBank = obj.funcBank;
            save(location,funcBank);
        end
        
        function [keys] = getKeys(obj)
            keys = fields(obj.funcBank);
        end
        
    end
    
    methods (Static)
        function [fileList] = generateFeatureFiles_forKey(fileList,oPath,key)
            for e = 1:numel(fileList)
                [pth nm ext] = fileparts(fileList{e}.fileName);
                fileList{e} = [oPath nm '--' key '.mat'];
            end
        end
        
        function [fileList] = generateFeatureFiles_forFile(file,oPath,keys)
            for e = 1:numel(keys)
                [pth nm ext] = fileparts(file.fileName);
                fileList{e} = [oPath nm '--' keys{e} '.mat'];
            end
        end
    end
end



