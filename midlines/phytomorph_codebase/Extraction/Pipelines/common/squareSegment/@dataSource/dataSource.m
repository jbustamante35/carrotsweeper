classdef dataSource < handle
    properties
        sourceList = {};
        loaderFunc = '';
        inMEM = 0;
        MEMstore = {};
    end
    
    methods
        
        function [obj] = dataSource(sourceList,loaderFunc,inMEM)
            obj.sourceList = sourceList;
            obj.loaderFunc = loaderFunc;
            obj.MEMstore = cell(numel(sourceList),1);
            obj.inMEM = inMEM;
        end
        
        % read from the stream
        function [tmpT] = read(obj,e)
            try
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % get the eth tensor
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if obj.inMEM && ~isempty(obj.MEMstore{e})
                    tmpT = obj.MEMstore{e};
                else
                    tmpT = obj.loaderFunc(obj.sourceList{e});
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % store the eth tensor
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if obj.inMEM
                    obj.MEMstore{e} = tmpT;
                end
            catch ME
                ME;
            end
        end
        
        % number of partices in stream
        function [z] = amount(obj)
            z = numel(obj.sourceList);
        end
        
        % pre-pressurize the dataSource
        function [] = fullLoad(obj)
            func = obj.loaderFunc;
            sourceList = obj.sourceList;
            for e = 1:numel(obj.sourceList)
                tmpT = func(sourceList{e});
                MEMstore{e} = tmpT;
                fprintf(['loaded:' num2str(e) ':' num2str(numel(obj.sourceList)) '\n'])
            end
            obj.MEMstore = MEMstore;
        end
    end
    
end

%{
    % to get perspectives on the dataSource one can apply a loader to
    obtain a draw, a facet,a filter
    loader: dataSource -n-> tensor
    facet: tensor --> vector(s)
    filter: vector(s) --> vector(s)


%}