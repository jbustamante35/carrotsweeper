classdef myV < handle
    %%% note not recursive complete !!!
    properties
        d; % data
    end    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        function [obj] = myV(varargin)
            if (nargin == 0)
                obj.d = [];
            else
                obj.d = varargin{1};
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set d
        function [obj] = setData(obj,d)
            obj.d = d;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get d
        function [d] = getData(obj)
            d = obj.d;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert to string
        function [] = toBSON(obj)
           obj.d = obj.getBSON();
        end        
        % convert to double
        function [] = toMATLAB(obj)
           obj.d = obj.getMATLAB();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert to BSON
        function [ret] = getBSON(obj)
            % if object is not in BSON state
            if ~strcmp(class(obj.d),'org.bson.BasicDBObject')
                import org.bson.*;
                import com.mongodb.*;
                %%%%%%%%%%%%%%
                %%% make output
                tensor = com.mongodb.BasicDBObject(numel(obj.d));
                % get the containerType
                containerType = class(obj.d);
                % switch on container type for index
                switch containerType
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'cell'
                        % use cell index
                        for e = 1:numel(obj.d)
                            % get the element
                            ele = obj.d{e};
                            classType = class(ele);
                            % switch on class type
                            switch classType
                                case 'char'
                                    tensor.put(num2str(e),java.lang.String(ele));
                                case 'single'
                                    tensor.put(num2str(e),java.lang.String(num2str(ele)));
                                case 'double'
                                    tensor.put(num2str(e),java.lang.String(num2str(ele)));
                                case 'uint8'
                                    tensor.put(num2str(e),java.lang.String(num2str(ele)));
                                otherwise
                                    try
                                        tensor.put(num2str(e),ele);
                                    catch ME
                                        fprintf(['error@myVar\n']);
                                    end
                            end
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'char'
                        % use array index
                        for e = 1:numel(obj.d)
                            tensor.put(num2str(e),java.lang.String(num2str(obj.d(e))));
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    otherwise % double,single,unit16,unit8 etc....
                        % use array index
                        for e = 1:numel(obj.d)
                            tensor.put(num2str(e),java.lang.String(num2str(obj.d(e))));
                        end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% make size
                sz = size(obj.d);
                tensor_sz =  com.mongodb.BasicDBObject();
                for e = 1:numel(sz)
                    tensor_sz.put(num2str(e),java.lang.String(num2str(sz(e))));
                end
                tensor_sz.put('ndims',java.lang.String(num2str(length(sz))));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make bson tensor
                ret = com.mongodb.BasicDBObject();
                ret.put('data',tensor);
                ret.put('size',tensor_sz);
                ret.put('container',containerType);                
            else
                ret = obj.d;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % convert to BSON - recursive
        function [ret] = getBSONr(obj)
            % if object is not in BSON state
            if ~strcmp(class(obj.d),'org.bson.BasicDBObject')
                import org.bson.*;
                import com.mongodb.*;
                %%%%%%%%%%%%%%
                %%% make output
                tensor = com.mongodb.BasicDBObject(numel(obj.d));
                % get the containerType
                containerType = class(obj.d);
                % switch on container type for index
                switch containerType
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'cell'
                        % use cell index
                        for e = 1:numel(obj.d)
                            % get the element
                            ele = obj.d{e};
                            classType = class(ele);
                            % switch on class type
                            switch classType
                                case 'char'
                                    tensor.put(num2str(e),java.lang.String(ele));
                                case 'single'
                                    tensor.put(num2str(e),java.lang.String(num2str(ele)));
                                case 'double'
                                    tensor.put(num2str(e),java.lang.String(num2str(ele)));
                                case 'uint8'
                                    tensor.put(num2str(e),java.lang.String(num2str(ele)));
                                otherwise
                                    try
                                        tensor.put(num2str(e),ele);
                                    catch ME
                                        fprintf(['error@myVar\n']);
                                    end
                            end
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'char'
                        % use array index
                        for e = 1:numel(obj.d)
                            tensor.put(num2str(e),java.lang.String(num2str(obj.d(e))));
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    otherwise % double,single,unit16,unit8 etc....
                        % use array index
                        for e = 1:numel(obj.d)
                            tensor.put(num2str(e),java.lang.String(num2str(obj.d(e))));
                        end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% make size
                sz = size(obj.d);
                tensor_sz =  com.mongodb.BasicDBObject();
                for e = 1:numel(sz)
                    tensor_sz.put(num2str(e),java.lang.String(num2str(sz(e))));
                end
                tensor_sz.put('ndims',java.lang.String(num2str(length(sz))));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make bson tensor
                ret = com.mongodb.BasicDBObject();
                ret.put('data',tensor);
                ret.put('size',tensor_sz);
                ret.put('container',containerType);                
            else
                ret = obj.d;
            end
        end
        % get the MATLAB
        function [ret] = getMATLAB(obj)
            if strcmp(class(obj.d),' com.mongodb.BasicDBObject')
                import org.bson.*;
                import com.mongodb.*;
                %%% get the data elements
                tensor = obj.d.get('data');
                tensor_sz = obj.d.get('size');
                containerType = obj.d.get('container');
                %%% make size of data
                dims = str2num(tensor_sz.get('ndims'));
                sz = zeros(1,dims);
                for d = 1:dims
                    sz(d) = str2double(tensor_sz.get(num2str(d)));
                end
                %%% preallocate data space - allocate at linear - reshape
                %%% later
                switch containerType
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'cell'
                        ret = cell(1,prod(sz));
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'char'
                        ret = '';
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    otherwise % double,unit8 etc....
                        ret = zeros(1,prod(sz));
                end
                %%% fill container
                switch containerType
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'cell'
                        % use cell index
                        for e = 1:prod(sz)
                            ret{e} = tensor.get(num2str(e));                            
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    case 'char'
                        % use array index
                        for e = 1:prod(sz)
                            ret(e) = tensor.get(num2str(e));
                        end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    otherwise % double,unit8 etc....
                        % use array index
                        for e = 1:prod(sz)
                            ret(e) = str2num(tensor.get(num2str(e)));
                        end
                end
                ret = reshape(ret,sz);
            else
                ret = obj.d;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%
        % recursive call to this
        function [out] = masterSwitch(in)
            switch class(in)
                case 'char'
                    out = convertChar(in);
                case 'Double'
                    out = convertTensor(in);
                case 'uint8'
                    out = convertTensor(in);
                case 'single'
                    out = convertTensor(in);
                case 'cell'
                    out = convertCell(in);
                case 'struct'
                    out = convertStruct(in)
            end
        end
        %%%%%%
        % convert cell
        function [out] = convertCell(in)
            out = BasicDBObject();            
            for e = 1:numel(in)
                out.append(fl,masterSwitch(in{e}));
            end            
        end        
        %%%%%%
        % convert struct
        function [out] = convertStruct(in)
            out = BasicDBObject();
            f = fieldnames(in);
            for fl = 1:numel(f)
                out.append(masterSwitch(in.f));
            end
        end
        %%%%%%
        % convert a matlab charater to String
        function [out] = convertChar(in)
            import java.lang.String.*;
            out = java.lang.String(in);
        end
        %%%%%%
        % convert tensor to 
        function [out] = convertTensor(in)
            import java.lang.Double.*;
            SZ = size(in);
            out = javaArray(java.lang.Double,prod(SZ));
            for e = 1:numel(in)
                out(e) = java.lang.Double(in(e));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
end