classdef tensorBundle < myHS_X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a set of tensor objects
    % this vs a set of sets or structs
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor  
        function [obj] = tensorBundle(varargin)
            % call super constructor - constrain to
            % tensor myT or lower
            obj = obj@myHS_X('myT');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert tensor in set nB as new trial
        function [] = insertB(obj,nB)
            for e = 1:numel(nB)
                obj{e}.insertT(nB{e});
            end
        end
        
    end
end
