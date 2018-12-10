classdef partialFunction
    properties
        func;
        algorithm;
    end
    
    methods
        function [obj] = partialFunction(func,varargin)
            if ~isa(func,'function_handle')
                tempdir = [tempname filesep];
                mkdir(tempdir);
                matFile = [tempdir func '.mat'];
                CMD = ['iget -f /iplant/home/nmiller/publicData/' func '.mat ' matFile];
                
                system(CMD);
                load(matFile);
            else
                obj.func = func;
                obj.algorithm = varargin{1};
            end
        end
        
        function [] = save(obj,filename)
            save(filename,'obj','-v7.3');
        end
        
        function [] = publish(obj)
            % render partial function to mat file under algorithm name
            tempdir = [tempname filesep];
            tempdir = ['/mnt/tetra/nate/' tempdir];
            mkdir(tempdir);
            matFile = [tempdir obj.algorithm '.mat'];
            obj.save(matFile);
            CMD = ['iput -f ' matFile ' /iplant/home/nmiller/publicData/'];
            system(CMD);
            
        end
        
        function [o] = evalFunction(obj,X)
            o = obj.func(X);
        end
        
    end
end

%{
    J = 3;
    g = @(X)X*J;
    f = partialFunction(g,'test');
    f.publish();
    h = partialFunction('test');
    l1 = f.evalFunction(4)
    l2 = h.evalFunction(4)
%}