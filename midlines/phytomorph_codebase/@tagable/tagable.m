classdef tagable < handle

    properties
        tags;
        transparent;
    end
    
    methods
        function [obj] = tagable()
            obj.tags = containers.Map();
        end
        
        
        function [] = setTag(obj,key,value,p)
            if nargin == 3
                p = true;
            end
            obj.tags(key) = value;
            if (p & obj.transparent)
                obj.renderTag(key);
            end
        end
        
        function [value] = getTag(obj,key)
            value = obj.tags(key);
        end
        
        function [] = removeTag(obj,key)
            obj.tags.remove(key);
            CMD = ['imeta rmw -d ' obj.remoteFileName ' ' key ' %'];
            [r,o] = system(CMD);
        end
        
        
        function [] = renderTags(obj)
            keys = obj.tags.keys;
            for k = 1:numel(keys)
                obj.renderTag(keys{k});
            end
        end
        
        function [] = renderTag(obj,key)
            CMD = ['imeta add -d ' obj.remoteFileName ' ' key ' ' obj.getTag(key)];
            if obj.exists()
                [r,o] = system(CMD);
            end
        end
        
        function [] = listTags(obj)
            keys = obj.tags.keys;
            for k = 1:numel(keys)
                v = obj.tags(keys{k});
                fprintf(['k:v-->' keys{k} ':' v '\n']);
            end
        end
        
    end
    
    
    methods (Abstract)
        
        iput(obj)
        
        iget(obj)
        
        existsR(obj)
        
        loadTags(obj)
        
    end
end