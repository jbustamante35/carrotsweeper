classdef CH
    properties
        % state information
        ch = [];
        % name information
        cn = [];
        % nummber of chromosome "units"
        cN = []
        % length of each unit
        cLv = [];
        % init type 
        type = '';
        % state range
        maxStates = [];
    end
    
    methods
        function [obj] = CH(type,maxStates,numberRange,lengthRange,numberDistribution,lengthDistribution)
            if nargin == 0
                type = 'rand';
                maxStates = 2;
                numberDistribution = 'Uniform';
                numberRange = {5 5};
                lengthDistribution = 'Uniform';
                lengthRange = {500 1000};
            end
            if nargin == 4
                numberDistribution = 'Uniform';
                lengthDistribution = 'Uniform';
            end
            para = CH.generateCHpara(numberDistribution,numberRange,lengthDistribution,lengthRange);
            obj.cN = para.cN;
            obj.cLv = para.cLv;
            [C] = generateCH(obj.cLv,maxStates,type);
            obj.ch = C.ch;
            obj.cn = C.cn;
        end
        
        
    end
    
    methods (Static)
        function [para] = generateCHpara(numberDistribution,numberRange,lengthDistribution,lengthRange)
            debug = 0;
            % generate the number of CH
            para.cN = generateCHn(numberDistribution,numberRange,{1});
            if debug;para.cN = 2;end
            % generate the length of each chr
            para.cLv = generateCHlength(lengthDistribution,lengthRange,para.cN);
            if debug;para.cLv = [100 100];end
        end
       
    end
end

%{

%}