classdef myiF < matlab.mixin.Copyable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % notes: added(13.01.25) - sven: there swhould be a discussion
    % reguarding the clustering of the p. in other words, the parameters
    % which are mapping from the domain to the codomain are themselves
    % a residing in a vector space. at some level, these vectors are siting
    % in the "functional" space (aka.hilbert space).  while these numbers 
    % may not full create a given functional "unit". these numbers will at
    % least, and maybe more, is a reduced version of the functional vector
    % in hilbert space. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%
    % "invert" functions
    properties
        f_handle;   % function handle
        if_handle;  % invert function handle
        p;          % parameters for forward function
        ip;         % parameters for invert function
    end
    
    %%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%
        % constructor
        function [obj] = myiF(varargin)
            obj.f_handle = [];
            obj.if_handle = [];
            if nargin >= 1
                obj.f_handle = varargin{1};
            end
            if nargin >= 2
                obj.if_handle = varargin{2};
            end
        end
        %%%%%%%%%%%%%
        % get parameters
        % notes: (13.01.25) the gP (get parameters) will create the
        % needed "numbers" to eval the function in the forward direction
        % each call to the gP will have access to the domain tensor and the
        % codomain tensor
        function [] = i(obj,d,c)
            obj.p = obj.if_handle(d,c,obj.ip);
        end
        %%%%%%%%%%%%%
        % eval
        function [ret] = e(obj,d)
            ret = obj.f_handle(d,obj.p);
        end
    end
    
    
end