classdef worker < matlab.mixin.Copyable
    
    properties
        jobQueue;
        kiosk;
    end
    
    methods
        
        function [obj] = worker()
            
        end
        
        function [] = attachJob(obj,jb)
            obj.jobQueue{end+1} = jb;
        end
        
        function [] = scanKioak(obj)
            
        end
        
        
        function [] = pwork(obj)
            for e = 1:numel(obj.jobQueue)
                try
                    obj.jobQueue{e}.run();
                catch
                    
                end
            end
            
        end
        
        function [] = work(obj)
            for e = 1:numel(obj.jobQueue)
                try
                    obj.jobQueue{e}.run();
                catch
                    
                end
            end
            
        end
    end
end