classdef phytoAflow_1 < myTb
    
    
    properties
    
    end
    
    methods
        
        function [obj] = phytoAflow_1(varargin)
            % super constructor - 1 rank fibre and 2 rank base
            obj = obj@myTb();
            obj.setAllRank(2,1);
            % set default views
            obj.view_props.props.LineStyle = '-';
            obj.view_props.props.Color = 'r';
            obj.view_props.type = 'phytoAsequence';
            % init point(s)
            if nargin == 1               
               obj.setData(varargin{1});
            end
        end
        
    end
    
end