classdef myProb < handle
    
    properties
        mu = [];
        sigma = [];
        
        X = [];
        ks = [];
        N = [];
        
        type = 'gaussian';
        gmm = [];
        ds = 1;
        nComp = 1;
    end
    
    methods
        function [obj] = myProb(mu,sigma)
            obj.mu = mu;
            obj.sigma = sigma;
        end
        
        function [] = fitToGMM(obj,data,nComp,ds)
            obj.type = 'gmm';
            obj.nComp = nComp;
            obj.ds = ds;
            obj.gmm = gmdistribution.fit(data(1:obj.ds:end,:),nComp);
            % removed for emergnce
            %obj.gmm = obj.gmm.preComp;
        end
        
        function [] = fitToKS(obj,data,X,N)
            obj.type = 'ksdensity';
            if nargin > 2
                obj.N = N;
                obj.ks = ksdensity(data,linspace(X(1),X(2),N));
                obj.X = linspace(X(1),X(2),N);
                obj.ks = obj.ks / sum(obj.ks);
            else
                
            end
        end
        
        function [sample] = drawSample(obj)
            if strcmp(obj.type,'gaussian');
                sample = mvnrnd(obj.mu,obj.sigma);
            else
                sample = random(obj.gmm,1);
            end
        end
        
        function [prob] = getProb(obj,observation,dn)
            if strcmp(obj.type,'gaussian')
                prob = dn*mvnpdf(observation',obj.mu,obj.sigma);
            elseif strcmp(obj.type,'ksdensity')
                prob = interp1(obj.X,obj.ks,observation);
            else
                prob = dn*pdf(obj.gmm,observation');
            end
        end
        
        function [] = update(obj,data)
            if strcmp(obj.type,'gaussian')
                obj.mu = mean(data,1);
                obj.sigma = cov(data);
            elseif strcmp(obj.type,'ksdensity')
                obj.fitToKS(data,obj.X,obj.N)
            else
                obj.gmm = gmdistribution.fit(data(1:obj.ds:end,:),obj.nComp);
                % removed for emergence
                %obj.gmm = obj.gmm.preComp;
            end
        end
    end
    
    
    
end