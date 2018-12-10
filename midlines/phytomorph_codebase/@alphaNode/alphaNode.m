classdef alphaNode < handle
    
    
    
    
    properties
        key;
        time;
        
        dataStore;
    end
    
    
    
    methods
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = alphaNode(dataStore,loadKey)
            % set datastore 
            obj.dataStore = dataStore;
            % if nargin
            [obj.key,obj.time] = alphaNode.generateKeyTime();
            
        end
        
        
        function [writePtr] = writeToDataStore()
            writePtr = obj.dataStore.getPtr();
            obj.dataStore.
        end
        
        function [] = removeFromDataStore()
            
        end
        
        function [] = loadFromDataStore()
        
        end
        
        
        
    end
    
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate key and time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [K,T] = generateKeyTime()
            % generate random index key
            K = randi((2^32)-1,1,'uint32');
            T = datenum(datetime);
        end
    end
    
end