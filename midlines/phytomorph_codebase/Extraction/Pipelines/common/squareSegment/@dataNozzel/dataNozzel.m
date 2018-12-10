classdef dataNozzle < handle
    
    properties
        dataSource;
        func;
        readsPerSource = 1;
        z = 0;
        U;
        C;
        tPtr = 1;
        sPtr = 1;
        ePtr = 1;
    end
    
    
    methods
        
        function [obj] = dataNozzel(func,dataSource,n)
            obj.dataSource = dataSource;
            obj.func = func;
            obj.readsPerSource = n;
            testRead = obj.read(1,1);
            obj.z = size(testRead,1);
        end
        
        
        
        function [d] = read(obj,e,n)
            d = obj.dataSource.load(e);
            d = obj.func(d,n);
        end
        
        function [b] = hasNext(obj)
            if obj.tPtr > obj.dataSource.size()*obj.readsPerSource
                b = 0;
            else
                b = 1;
            end
        end
        
        function [] = resetPointer(obj)
            obj.ePtr = 1;
            obj.sPtr = 1;
            obj.tPtr = 1;
        end
        
        function [d] = next(obj)
            d = read(obj,obj.ePtr,obj.sPtr);
            obj.ePtr = obj.ePtr + 1;
            obj.sPtr = obj.sPtr + 1;
            if obj.sPtr > obj.readsPerSource
                obj.sPtr = 1;
            end
            obj.tPtr = obj.tPtr + 1;
        end
        
        function [U] = expectedValue(obj)
            if isempty(obj.U)
                s = zeros(obj.z,1);
                z = [];
                while obj.hasNext()
                    I = obj.next();
                    s = s + sum(I,2);
                    z = [z size(I,2)];
                end
                obj.resetPointer();
                % get expected value
                obj.U = s * sum(z)^-1;
            end
            U = obj.U;
        end
        
        function [C] = covMatrix(obj)
            if isempty(obj.C)
                obj.C = zeros(obj.z,obj.z);
                U = obj.expectedValue();
                z = 0;
                while obj.hasNext()
                    % get next from stream
                    I = obj.next();
                    % subtract mean
                    I = bsxfun(@minus,I,U);
                    % create COV
                    obj.C = obj.C + I*I';
                    % new image size
                    z = z + size(I,2);
                end
                obj.resetPointer();
                obj.C = obj.C * z.^-1;
            end
            C = obj.C;
        end
        
        
        
    end
end