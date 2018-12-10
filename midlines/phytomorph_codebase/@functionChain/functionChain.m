classdef functionChain
    properties
    end
    
    methods
        function [obj] = functionChain(dbFile)
            % if no file name is given-then default to in memory database
            if nargin == 0
                dbFile = ':memory:';
            end
            % open database
            mksqlite(0,'open',dbFile);
            % set blob type
            mksqlite('typedBLOBs', 2 );
            % wrapping
            mksqlite('param_wrapping', 0 );
            % create table for args
            mksqlite('CREATE TABLE fArgs (args)');
            % create table for function strings
            mksqlite('CREATE TABLE functions (func)');
        end
        
        
        function [] = putFunction(obj,fIndex,func,args)
            mksqlite('INSERT INTO fArgs (args) VALUES (?)', args);
            func = func2str(func)
            mksqlite('INSERT INTO functions (func) VALUES (?)', func);
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
    end
    
    
end

%{
    C = functionChain()
    func = @(X,Y)plus(X,Y);
    args = {4};
    C.putFunction(1,func,args)
    out = C.apply(1)


    C = functionChain();
    func = @(X,alongDim,numberFreq)extractionAtom0(X,alongDim,numberFreq);
    args = {1,20};
    C.putFunction(1,func,args);
    out = C.apply(rand(100,100));

%}