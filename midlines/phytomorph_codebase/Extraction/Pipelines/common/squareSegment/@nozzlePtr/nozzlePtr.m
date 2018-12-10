classdef nozzlePtr < handle
    properties
        currentValue = 1;
        majorValue = 1;
        minorValue = 1;
        majorMax;
        minorMax;
    end
    
    methods
        function [obj] = nozzlePtr(majorMax,minorMax)
            obj.majorMax = majorMax;
            obj.minorMax = minorMax;
        end
        
        function [] = reset(obj)
            obj.majorValue = 1;
            obj.minorValue = 1;
            obj.currentValue = 1;
        end
        
        function [majorValue,minorValue] = next(obj)
            
            
            
            [obj.majorValue,obj.minorValue] = obj.calculateNindex(obj.currentValue);
            obj.currentValue = obj.currentValue + 1;
            majorValue = obj.majorValue;
            minorValue = obj.minorValue;
            %{
            
            obj.minorValue = obj.minorValue + 1;
            if obj.minorValue > obj.minorMax
                obj.minorValue = 1;
            end
            obj.majorValue = obj.majorValue + 1;
            %}
            
            
            
        end
        
        function [majorValue,minorValue] = calculateNindex(obj,value)
            majorValue = floor(value/obj.minorMax);
            if majorValue == 0
                majorValue = 1;
            end
            minorValue = mod(value-1,obj.minorMax)+1;
        end
        
        function [ret] = hasNext(obj)
            if obj.currentValue > obj.majorMax*obj.minorMax
                ret = 0;
            else
                ret = 1;
            end
        end
        
        
    end
end