classdef nozzleManifold < handle
    properties
        nozzlePtr;
        sourceNozzles;
        manifoldFunc;
        readsPerSource;
        z;
    end
    
    methods
        function [obj] = nozzleManifold(sourceNozzles,manifoldFunc,readsPerSource)
            obj.sourceNozzles = sourceNozzles;
            obj.manifoldFunc = manifoldFunc;
            obj.readsPerSource = readsPerSource;
            obj.nozzlePtr = nozzlePtr(obj.amount()/readsPerSource,readsPerSource);
            % peak read for object size
            testRead = obj.read(1,1);
            obj.z = size(testRead,1);
        end
        
        function [d] = next(obj)
            [majorValue,minorValue] = obj.nozzlePtr.next();
            %d = read(obj,majorValue,minorValue);
            d = obj.read(majorValue,minorValue);
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
        
         % nessary from parent class
        function [d] = read(obj,varargin)
            
            if numel(varargin) == 1
                % fraction the request into major and minor
                [majorValue,minorValue] = obj.nozzlePtr.calculateNindex(varargin{1});
            else
                majorValue = varargin{1};
                minorValue = varargin{2};
            end
            
            for s = 1:numel(obj.sourceNozzles)
                % read the eth object from the stream
                d{s} = obj.sourceNozzles{s}.read(majorValue);
            end
            
            
            % fraction the eth object 
            d = obj.manifoldFunc(d,majorValue,minorValue);
        end
        
        function [ret] = hasNext(obj)
            ret = obj.nozzlePtr.hasNext();
            %{
            for s = 1:numel(obj.sourceNozzles)
                ret(s) = obj.sourceNozzles{s}.hasNext();
            end
            ret = all(ret);
            %}
        end
        
        function [] = resetPtr(obj)
            obj.nozzlePtr.reset();
            %{
            for s = 1:numel(obj.sourceNozzles)
                obj.sourceNozzles{s}.resetPtr();
            end
            %}
        end
        
        function [a] = amount(obj)
            for s = 1:numel(obj.sourceNozzles)
                a(s) = obj.sourceNozzles{s}.amount();
            end
            a = mean(a)*obj.readsPerSource;
        end
        
        function [C] = collectFullSpray(obj)
            
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
            
            %{
            C = [];
            cnt = 1;
            totalParticles = obj.amount();
            while obj.hasNext()
                tic;
                fprintf(['start reading particle(s):' num2str(cnt) ':' num2str(totalParticles) '\n']);
                C = [C obj.next()];
                fprintf(['end reading particle(s):' num2str(cnt) ':' num2str(totalParticles) ':' num2str(toc) '\n']);
                cnt = cnt + 1;
            end
            obj.resetPtr();
            %}
            
        end
    end
end