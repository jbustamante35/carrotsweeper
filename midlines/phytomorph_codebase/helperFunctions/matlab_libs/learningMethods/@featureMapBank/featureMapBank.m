classdef featureMapBank < handle
    properties
        oPath = '';
        featureFunctionBank;
        typeBank;
        % bank of features listed by type : 
        %       1) hash type : one to many
        %       2) moniker type : one to many
        bank;
    end
    
    methods
        function [obj] = featureMapBank(varargin)
            if nargin == 1
                
            end
        end
        
        function [] = depositFeatureMap(obj,featureMapObject)
            typeKey = featureMapObject.getTypeKey();
            obj.bank.byType.objectList.(typeKey){end+1} = featureMapObject;
            depositKey = ['k_' num2str(string2hash(featureMapObject.key))];
            obj.bank.byID.objectList.(depositKey) = featureMapObject;
        end
        
        function [rfmo] = withdrawFeatureMap(obj,fmoKey)
            withDrawKey = ['k_' num2str(string2hash(fmoKey))];
            rfmo = obj.bank.byID.objectList.(withDrawKey);
        end
        
        function [] = registerType(obj,fmo)
            if ~isfield(obj.typeBank,fmo.typeKey)
                obj.typeBank.(fmo.typeKey) = fmo;
                obj.bank.byType.objectList.(fmo.typeKey) = {};
                obj.bank.byType.monikerList.(fmo.typeKey) = fmo.moniker;
            end
        end
        
        function [monikerList] = getFeaturesByType(obj)
            f = fields(obj.bank.byType.monikerList);
            for e = 1:numel(f)
                monikerList{e} = obj.bank.byType.monikerList.(f{e});
            end
        end
        
        function [featureMapObjects] = getFeatureMap(obj,search_moniker)
            featureMapObjects = {};
            f = fields(obj.bank.byType.monikerList);
            for e = 1:numel(f)
                monikerToSearch = obj.bank.byType.monikerList.(f{e});
                if ~isempty(findstr(search_moniker,monikerToSearch))
                    featureMapObjects = obj.bank.byType.objectList.(f{e});
                end
            end
        end
    end
    
    methods (Static)
        function [index] = inverseSelect(selectFunction,classifyFeatures,transformFeatures,featureMapObjects)
            for e = 1:numel(featureMapObjects)
                tmpD = featureMapObjects{e}.getData();
                tmpD = transformFeatures(tmpD);
                tmpD = classifyFeatures(tmpD);
                index{e} = selectFunction(tmpD);
            end
        end
        
        function [indexSets] = intersectIndex(indexSets,featureMapObjectSets)
            for e = 1:numel(indexSets)
                for m = 1:numel(featureMapObjectSets)
                    indexSets{e} = intersect(indexSets{e},featureMapObjectSets{m}{e}.getIndexMap());
                end
            end
        end
       
        function [data] = loadFeatureMapsAtAcrossIndexSets(featureMapObjectSets,indexSets)
            data = [];
            % for each request "type" request
            if nargin == 1
                 [ptToArray] = featureMapBank.calcAllocate(featureMapObjectSets);
                 fprintf(['Start preallocate feature map load \n']);tic
                 data = zeros(ptToArray{end,end}(2),ptToArray{end,end}(4));
                 fprintf(['Done preallocate feature map load:' num2str(toc) '\n']);
            end
            for e = 1:numel(featureMapObjectSets)
                tmpD = [];
                % for each map in the "type"
                for m = 1:numel(featureMapObjectSets{e})
                    tm = clock;
                    if nargin > 1
                        if iscell(indexSets)
                            tmpD = cat(2,tmpD,featureMapObjectSets{e}{m}.getData(indexSets{m}));
                        else
                            tmpD = cat(2,tmpD,featureMapObjectSets{e}{m}.getData(indexSets));
                        end
                    else
                        % speed up via pre-allocation of memory in the case
                        % of total feature maps reads
                        str1 = ptToArray{e,m}(1);
                        stp1 = ptToArray{e,m}(2);
                        str2 = ptToArray{e,m}(3);
                        stp2 = ptToArray{e,m}(4);
                        fprintf(['Start disk read feature map load \n']);tic
                        tmpD = featureMapObjectSets{e}{m}.getData();
                        fprintf(['Done disk read feature map load:' num2str(toc) '\n']);
                        fprintf(['Start store tmp read feature map load \n']);tic
                        data(str1:stp1,str2:stp2) = tmpD;
                        fprintf(['Done store tmp read feature map load:' num2str(toc) '\n']);
                    end
                    fprintf(['Done loading map: ' num2str(m) ':feature:' num2str(e) '@' num2str(etime(clock,tm)) '\n']);
                end
                if nargin ~= 1
                    data = cat(1,data,tmpD);
                end
            end
        end
        
        function [ptToArray] = calcAllocate(featureMapObjectSets)
            % for each request "type" request
            str1 = 1;
            for e = 1:numel(featureMapObjectSets)
                % for each map in the "type"
                str2 = 1;
                for m = 1:numel(featureMapObjectSets{e})
                    sz = featureMapObjectSets{e}{m}.getFeatureMapSize();
                    ptToArray{e,m} = [str1 str1 + sz(1) - 1 str2 str2 + sz(2) - 1];
                    str2 = str2 + sz(2);
                end
                str1 = str1 + sz(1);
            end
        end
    end
    
end