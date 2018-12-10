classdef FRE < handle
    
    properties
        f_func;
        r_func;
        e_func;
    end
    
    
    methods
        
        function [obj] = FRE()
        
        end
        
        function [] = setF(obj,f_func)
            obj.f_func = f_func;
        end
        
        function [] = setR(obj,r_func)
            obj.r_func = r_func;
        end
        
        function [] = setE(obj,e_func)
            obj.e_func = e_func;
        end
        
        function [] = run(obj,in,subIndex)
            n = obj.f_func.getNargout();
            
            for e = 1:numel(in)
                [out{1:n}] = obj.f_func.run(in{e});
                %f{e} = obj.f_func.run(in{e});
            end
            
            
            [data] = featureMapBank.loadFeatureMapsAtAcrossIndexSets({f},subIndex);
            bi = obj.r_func(data);
            
            
            
            for e = 1:numel(f)
                er{e} = obj.e_func.run(f{e},bi);
            end
            
        end
    end
end