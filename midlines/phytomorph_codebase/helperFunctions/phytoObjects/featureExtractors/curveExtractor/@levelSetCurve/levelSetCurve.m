classdef levelSetCurve  < curveExtractor
    
    properties
      
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = levelSetCurve(varargin)
            % super constructor
            obj = obj@curveExtractor();
            % hookin level set extractor
            obj.funcList{1} = @(X,l)contourc(X,l);
            % add number of level sets
            if nargin == 1
                obj.para = varargin{1};
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract curves
        function [cb] = extractCurves(obj,X)
            cb = myHS_X('phytoAcurve');
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % reader
            X = myReader(X,'toGray');
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % run func
            c = obj.funcList{1}(X,obj.para);   
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % init vars
            str = 1;
            len = c(2,str);
            en = str + len;
            cnt = 1;
            stop = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop over curves
            while ~stop
                % get current curve
                cb{cnt} = phytoAcurve(c(:,str+1:en));
                data = c(:,str+1:en);
                cb{cnt} = phytoAcurve(data');
                % increment counter
                cnt = cnt + 1;
                % if end is at stop then flag
                if en == size(c,2)
                    stop = 1;
                else
                    str = en + 1;
                    len = c(2,str);
                    en = str + len;    
                end
            end
            
            
        end
    end
end
