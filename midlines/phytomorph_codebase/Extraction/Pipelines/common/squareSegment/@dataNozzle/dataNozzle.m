classdef dataNozzle < handle
    
    properties
        dataSource;
        func;
        
        nozzlePtr;
        
        
        readsPerSource = 1;
        
        
        z = 0;
        U;
        C;
        
        
    end
    
    
    methods
        
        function [obj] = dataNozzle(func,dataSource,n)
            try
                obj.dataSource = dataSource;
                obj.func = func;
                obj.readsPerSource = n;
                obj.nozzlePtr = nozzlePtr(dataSource.amount(),n);

                % peak read for object size
                testRead = obj.read(1,1);
                obj.z = size(testRead,1);
            catch ME
                getReport(ME)
            end
        end
        
        
       
        % number of objects available from the source
        function [z] = amount(obj)
            z = obj.dataSource.amount();
            z = z*obj.readsPerSource;
        end
        
        % dimensions
        
        
        
        % nessary from parent class
        function [d] = read(obj,varargin)
            fprintf(['calling read on: ' func2str(obj.func) '\n'])
            if numel(varargin) == 1
                % fraction the request into major and minor
                [majorValue,minorValue] = obj.nozzlePtr.calculateNindex(varargin{1});
            else
                majorValue = varargin{1};
                minorValue = varargin{2};
            end
            
            % read the eth object from the stream
            d = obj.dataSource.read(majorValue);
            % fraction the eth object 
            d = obj.func(d,majorValue,minorValue);
        end
        
        
        
        function [b] = hasNext(obj)
            b = obj.nozzlePtr.hasNext();
        end
        
        function [] = resetPtr(obj)
            obj.nozzlePtr.reset();
        end
        
        function [d] = next(obj)
            [majorValue,minorValue] = obj.nozzlePtr.next();
            d = read(obj,majorValue,minorValue);
        end
        
        function [C]= collectFullSpray(obj)
            C = [];
            cnt = 1;
            totalParticles = obj.amount();
            countGauge = nozzleGauge(obj,{@(X,v)particleCountingGauge(X,v)});
            particleCount = countGauge.takeMetric();
            C = zeros(obj.z,particleCount{1});
            ptr = 1;
            while obj.hasNext()
                tic;
                fprintf(['start reading particle(s):' num2str(cnt) ':' num2str(totalParticles) '\n']);
                tmp = obj.next();
                stp = ptr + size(tmp,2)-1;
                C(:,ptr:stp) = tmp;
                ptr = stp + 1;
                %C = [C obj.next()];
                fprintf(['end reading particle(s):' num2str(cnt) ':' num2str(totalParticles) ':' num2str(toc) '\n']);
                cnt = cnt + 1;
            end
            obj.nozzlePtr.reset();
        end
        
        
        
        function [tNozzle] = getTransductionNozzle(obj,transductionFunction)
            % collect the particles
            C = obj.collectFullSpray();
            % evaluate the transduction (fitting) functor
            func = transductionFunction(C);
            % construct func for loading into new nozzle
            %func = @(X,e0,e1)predict(evaulatedTransducor,X');
            % construct and load nozzle
            tNozzle = dataNozzle(func,obj,1);
        end
        
        
        function [rNozzle] = getReductionNozzle(obj,k)
            U = obj.expectedValue();
            C = obj.covMatrix();
            [E,D] = eigs(C,k);
            func = @(X,e0,e1)reductionNozzle(X,U,E);
            rNozzle = dataNozzle(func,obj,1);
        end
        
        % gauages on the stream
        function [U] = expectedValue(obj)
            % obtaining expected value
            fprintf(['starting calc on expected value \n']);
            if isempty(obj.U)
                s = zeros(obj.z,1);
                z = [];
                totalParticles = obj.amount();
                cnt = 1;
                while obj.hasNext()
                    tic
                    fprintf(['start reading particle:' num2str(cnt) ':' num2str(totalParticles) '\n']);
                    I = obj.next();
                    s = s + sum(I,2);
                    z = [z size(I,2)];
                    fprintf(['end reading particle:' num2str(cnt) ':' num2str(totalParticles) ':' num2str(toc) '\n']);
                    cnt = cnt + 1;
                end
                obj.nozzlePtr.reset();
                % get expected value
                obj.U = s * sum(z)^-1;
            end
            U = obj.U;
            fprintf(['ending calc on expected value \n']);
        end
        
        function [C] = covMatrix(obj)
            if isempty(obj.C)
                obj.C = zeros(obj.z,obj.z);
                U = obj.expectedValue();
                z = 0;
                totalParticles = obj.amount();
                cnt = 1;
                while obj.hasNext()
                    tic
                    fprintf(['start reading particle:' num2str(cnt) ':' num2str(totalParticles) '\n']);
                    % get next from stream
                    I = obj.next();
                    % subtract mean
                    I = bsxfun(@minus,I,U);
                    % create COV
                    obj.C = obj.C + I*I';
                    % new image size
                    z = z + size(I,2);
                    fprintf(['end reading particle:' num2str(cnt) ':' num2str(totalParticles) ':' num2str(toc) '\n']);
                    cnt = cnt + 1;
                end
                 obj.nozzlePtr.reset();
                obj.C = obj.C * z.^-1;
            end
            C = obj.C;
        end
        
        
        
    end
end