classdef ptde <  handle
    
    properties
        oPath = [];
        name = '';
        featureMapBank_file = '';
        featureMapBank;
        featureFunctionBank_file = '';
        featureFunctionBank;
    end
    
    methods
    
        function [obj] = ptde(oPath,name)
            % set the path for the point detector
            obj.oPath = oPath;
            obj.name = name;
            this_fileName = [obj.oPath obj.name '.mat'];
            if exist(this_fileName)
                load(this_fileName,'obj');
            else
                %%% feature map bank
                % set the file for the feature map bank
                obj.featureMapBank_file = [oPath 'featureMapBank.mat'];
                if exist(obj.featureMapBank_file)

                else
                    obj.featureMapBank = featureMapBank();
                end
                %%% feature function bank
                % set the file for the feature map bank
                obj.featureFunctionBank_file = [oPath 'featureFunctionBank.mat'];
                if exist(obj.featureFunctionBank_file)

                else
                    obj.featureFunctionBank = featureFunctionBank();
                end
            end
        end
        
        % add a function to the function bank
        function [f_key] = addFunction(obj,funcObject)
            f_key = obj.featureFunctionBank.addFunction(funcObject);
        end
        
        function [] = depositFeatureMaps(obj,featureMapObject)
            obj.featureMapBank.depositFeatureMap(featureMapObject);
        end
        
        function [rfmo] = withdrawFeatureMap(obj,fmoKey)
            rfmo = obj.featureMapBank.withdrawFeatureMap(fmoKey);
        end
        
        function [] = registerFeatureMapType(obj,featureMap)
            obj.featureMapBank.registerType(featureMap);
        end
        
        function [] = persist(obj)
            oFile = [obj.oPath obj.name '.mat'];
            save(oFile,'obj');
        end
        
        function [monikerList] = getFeatureMapMonikerList(obj)
            [monikerList] = obj.featureMapBank.getFeaturesByType();
        end
        
        function [featureMapObjects] = getFeatureMaps_byMonikerSearch(obj,search_moniker)
            [featureMapObjects] = obj.featureMapBank.getFeatureMap(search_moniker);
        end
        
      
    end
end