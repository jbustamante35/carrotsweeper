classdef logisticFit < modelFit
    
    properties
    
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = logisticFit()
        
        end
        %%%%%%%%%%%%%%%%%%%
        % fit the data
        function [] = fit(obj,X,Y)
            if nargin == 1
                X = obj.Xdata;
                Y = obj.Xdata;
            else
                obj.Xdata = X;
                obj.Ydata = Y;
            end
            
            warning off
            % init the sweep values
            P{1} = linspace(.5,8,3);
            P{2} = linspace(.001,.01,3);        
            % init the values
            p1 = P{1};
            p2 = P{2};
            cnt = 1;
            [J sidx] = sort(Y);
            vmax = mean(Y(end-2));
            for i = 1:size(p1,2)
                for j = 1:size(p2,2)
                    % render a set of values for defaults
                    p0 = [mean(Y),vmax,p1(i),p2(i),mean(X)]; % init Parameters             
                    % set the fixed and suggested
                    suggested.k = ones(size(p0));
                    suggested.v = zeros(size(p0));
                    fixed.k = ones(size(p0));
                    fixed.v = zeros(size(p0));
                    % infuse suggestions
                    p0 = p0.*suggested.k + suggested.v;
                    options = statset('Robust','on','WgtFun','logistic');
                    try
                        [p0,J1,J2,J3,mse(cnt)] = nlinfit(X,Y,@(X,Y)veloc_SPEC(X,Y,fixed),p0,options);
                        p0 = p0.*fixed.k + fixed.v;
                        beta(cnt,:) = p0;
                        cnt = cnt + 1;
                    catch ME
                        cnt = cnt + 1;
                        beta(cnt,:) = p0;
                        mse(cnt) = inf;
                    end
                end
            end
            % find complex numbers for error and set to inf
            fidx = find(imag(mse) ~= 0);
            fidx1 = find(all(imag(beta) ~= 0,2));
            mse(fidx) = inf;
            mse(fidx1) = inf;
            % select min error
            [JUNK,sidx] = min(mse);
            obj.parameters = beta(sidx,:);
            
        end
        %%%%%%%%%%%%%%%%%%%
        % eval the logistic function
        function [Y] = fnval(obj,X)
            Y = veloc_SPEC(obj.parameters,X);
        end
    end
    
end