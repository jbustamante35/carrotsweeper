classdef xM < handle
    properties
        % cross over matrix
        m = [];
        % number of cross over sites
        xP = [];
        % 
        xSites = {};
        xVec = {};
        exeProb = {};
        cLv = [];
        crossOverSitesNumberRange = {};
        crossOverSitesNumberDistribution = '';
        crossOverSitesDistribution = '';
    end
    methods
        function [obj] = xM(chr,crossOverSitesNumberRange,crossOverSitesNumberDistribution,crossOverSitesDistribution)
            if nargin == 1
                crossOverSitesNumberRange = {2 10};
                crossOverSitesNumberDistribution = 'Uniform';
                crossOverSitesDistribution = 'Uniform';
            end
            obj.crossOverSitesNumberRange = crossOverSitesNumberRange;
            obj.crossOverSitesNumberDistribution = crossOverSitesNumberDistribution;
            obj.crossOverSitesDistribution = crossOverSitesDistribution;
            %%%%%%%%%%%%%%%%%%%%%%
            % generate  crossover matrix
            for e = 1:numel(chr.cLv)
                debug = 0;
                %%%%%%%%%%%%%%%%%%%%%%
                % generate number of xOver sites
                obj.xP(e) = generateNumberOfXoverSites(crossOverSitesNumberDistribution,crossOverSitesNumberRange);
                if debug;obj.xP(e) = 2;end
                %%%%%%%%%%%%%%%%%%%%%%
                % generate the xover sites
                [obj.xSites{e} obj.xVec{e}] = generateXoverLocs(crossOverSitesDistribution,{1 obj.xP(e)},{1 chr.cLv(e)});
                %%%%%%%%%%%%%%%%%%%%%%
                % generate the threshold for xovr to happen
                [obj.exeProb{e}] = generateExecuteProb(numel(obj.xSites{e}));
                if debug;obj.exeProb{e} = zeros(size(obj.exeProb{e}));end
            end
            obj.cLv = chr.cLv;
        end
        
        function [] = gen_xM(obj)
            obj.m = generateXover(obj.xSites,obj.cLv,obj.exeProb);
        end
    end
end