classdef transformationList < globalDB
    properties
    end
    
    methods
        function [obj] = transformationList(varargin)
            obj = obj@globalDB(varargin{:});
            % create table for args
            mksqlite(obj.id,'CREATE TABLE transformationList (name,tforms)');
        end
        
        
        function [] = putTransformation(obj,name,tform)
            mksqlite(obj.id,'INSERT INTO transformationList VALUES (?,?)', {name,tform});
        end
        
        
        function [out] = apply(obj,X,stopN)
            out = 1;
            args = mksqlite( 'SELECT args FROM fArgs' );
            funcs = mksqlite( 'SELECT func FROM functions' );
            if nargin == 2;stopN = numel(args);end
            for e = 1:min(stopN,numel(args))
                tmpFunc = str2func(funcs(e).func);
                X = tmpFunc(X,args(e).args{:});
            end
        end
        
        function [] = close(obj)
            mksqlite(obj.id,'close')
        end
    end
    
    
end
%{

    C = transformationList();
    tF = rand(100,2);
    C.putTransformation('pixelLocations',tF);

    func = @(X,Y)plus(X,Y);
    args = {4};
    C.putFunction(1,func,args)
    out = C.apply(1)
%}