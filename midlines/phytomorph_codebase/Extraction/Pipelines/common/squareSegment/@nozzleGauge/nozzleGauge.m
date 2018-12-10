classdef nozzleGauge < handle
    properties
        gaugeFunc;
        nozzle;
        M;
    end
    
    methods
        
        
        
        function [obj] = nozzleGauge(nozzle,gaugeFuncs)
            obj.gaugeFunc = gaugeFuncs;
            obj.nozzle = nozzle;
            obj.M = cell(numel(gaugeFuncs),1);
        end
        
        
        
        function [M] = takeMetric(obj)
            try
                for m = 1:numel(obj.gaugeFunc)
                    EM(m) = isempty(obj.M{m});
                end

                if any(EM)
                    M = cell(size(obj.M));
                    while obj.nozzle.hasNext()
                        tmpParticles = obj.nozzle.next();
                        for m = 1:numel(obj.gaugeFunc)
                            M{m} = obj.gaugeFunc{m}(tmpParticles,M{m});
                        end
                    end
                    obj.nozzle.resetPtr();
                    obj.M = M;
                end
                M = obj.M;
            catch ME
                ME;
            end
        end
        
        
        
        
    end
end