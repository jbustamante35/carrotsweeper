function [X,map] = readtif(filename, varargin)
%READTIF Read an image from a TIFF file.
%   [X,MAP] = READTIF(FILENAME) reads the first image from the
%   TIFF file specified by the string variable FILENAME.  X will
%   be a 2-D uint8 array if the specified data set contains an
%   8-bit image.  It will be an M-by-N-by-3 uint8 array if the
%   specified data set contains a 24-bit image.  MAP contains the
%   colormap if present; otherwise it is empty. 
%
%   [X,MAP] = READTIF(FILENAME, N, ...) reads the Nth image from the
%   file.
%
%   READTIF accepts trailing parameter-value pairs.  The parameter names and
%   corresponding values are:
%
%       Parameter Name        Value
%       --------------        -----
%       Index                 Scalar; same as input parameter N above
%
%       Info                  Struct; same as input parameter INFO above
%
%       PixelRegion           {ROWS, COLS}
%                             reads a region of pixels from the file. ROWS
%                             and COLS are two- or three-element vectors,
%                             where the first value is the start location,
%                             and the last value is the ending location. In
%                             the three-value syntax, the second value is the
%                             increment.
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Copyright 1984-2011 The MathWorks, Inc. 
%   $Revision: 1.1.6.14 $  $Date: 2011/05/17 02:28:09 $

args = parse_args(varargin{:});
checkIndex(args.index)
args.filename = filename;
args.pixelregion = process_region(args.pixelregion);

if ~isempty(args.info) && isfield(args.info, 'Offset')
    % If INFO input is provided, use the Offset field to determine where the
    % specified image is located in the file.
    
    if args.index > numel(args.info)
        error(message('MATLAB:imagesci:readtif:indexInfoMismatch'));
    end

    args.offset = args.info(args.index).Offset;
end

[X, map, details] = rtifc(args);

map = double(map)/65535;

if (details.Photo == 8)
    % TIFF image is in CIELAB format.  Issue a warning that we're
    % converting the data to ICCLAB format, and correct the a* and b*
    % values.
    
    % First, though, check to make sure we have the expected number of
    % samples per pixel.
    if (size(X,3) ~= 1) && (size(X,3) ~= 3)
        error(message('MATLAB:imagesci:readtif:unexpectedCIELabSamplesPerPixel'));
    end
    
    % Now check that we have uint8 or uint16 data.
    if ~ (isa(X,'uint8') || isa(X,'uint16'))
        error(message('MATLAB:imagesci:readtif:wrongCieLabDatatype', class( X )));
    end
    
    warning(message('MATLAB:imagesci:readtif:CielabConversion'));
    
    X = cielab2icclab(X);
elseif (details.Photo == 0)
    if (details.BitsPerSample == 1)
        % TIFF image is MinIsWhite and logical.  Invert it.
        X = ~X;
    elseif (details.Photo == 0) && (details.BitsPerSample <= 8)
        % TIFF image is MinIsWhite and not more than one byte.  Invert it.
        X = 2^(details.BitsPerSample) - 1 - X;
    end
end



function args = parse_args(varargin)
%PARSE_ARGS  Convert input arguments to structure of arguments.

args.index = 1;
args.pixelregion = [];
args.info = [];

params = {'index', 'pixelregion', 'info','workingDir'};

p = 1;
while (p <= nargin)
    
    argp = varargin{p};
    if (isnumeric(argp))

        args.index = argp;
        p = p + 1;
        
    elseif (ischar(argp))
        
        idx = find(strncmpi(argp, params, numel(argp)));
        
        if (isempty(idx))
            error(message('MATLAB:imagesci:readtif:unknownParam', argp))
        elseif (numel(idx) > 1)
            error(message('MATLAB:imagesci:readtif:ambiguousParam', argp))
        end
        
        if (p == nargin)
            error(message('MATLAB:imagesci:readtif:missingValue', argp))
        end
        
        args.(params{idx}) = varargin{p + 1};
        p = p + 2;
        
    else
        
        error(message('MATLAB:imagesci:readtif:paramType'))
        
    end
            
end

check_info(args.info);

function check_info(info)
%CHECK_INFO Issue error message if user passed in invalid info struct.

if isempty(info)
    return
end

if ~all(isfield(info, {'Filename', 'FileModDate', 'FileSize', ...
                       'Format', 'FormatVersion', 'Width', ...
                       'Height', 'BitDepth', 'ColorType', ...
                       'FormatSignature'}))
    error(message('MATLAB:imagesci:readtif:invalidInfoStruct'));
end


function region_struct = process_region(region_cell)
%PROCESS_PIXELREGION  Convert a cells of pixel region info to a struct.

region_struct = struct([]);
if isempty(region_cell)
    % Not specified in call to readtif.
    return;
end

if ((~iscell(region_cell)) || (numel(region_cell) ~= 2))
    error(message('MATLAB:imagesci:readtif:pixelRegionCell'))
end

for p = 1:numel(region_cell)
    
    checkIntegers(region_cell{p});
    
    if (numel(region_cell{p}) == 2)
        
        start = max(0, region_cell{p}(1) - 1);
        incr = 1;
        stop = region_cell{p}(2) - 1;
        
    elseif (numel(region_cell{p}) == 3)
        
        start = max(0, region_cell{p}(1) - 1);
        
        if (~isinf(region_cell{p}(2)))
            incr = region_cell{p}(2);
        else
            error(message('MATLAB:imagesci:readtif:infIncrement'));
        end
        
        stop = region_cell{p}(3) - 1;
       
    else
        
        error(message('MATLAB:imagesci:readtif:tooManyPixelRegionParts'));
        
    end
        
    if (start > stop)
        error(message('MATLAB:imagesci:readtif:badPixelRegionStartStop'))
    end

    if (incr < 1)
        error(message('MATLAB:imagesci:readtif:badPixelRegionIncrement'))
    end

    region_struct(p).start = start;
    region_struct(p).incr = incr;
    region_struct(p).stop = stop;

end


function checkIntegers(inputValue)

if (~isnumeric(inputValue) || ...
    any((rem(inputValue, 1) ~= 0) & ~(isinf(inputValue))) || ...
    any(inputValue < 0))
    
    error(message('MATLAB:imagesci:readtif:regionPartNotNumeric'));

end


function checkIndex(inputIndex)

if (~isPositiveFiniteIntegerScalar(inputIndex))
    
    error(message('MATLAB:imagesci:readtif:badIndex'))

end


function tf = isPositiveFiniteIntegerScalar(values)

tf = isnumeric(values) && ...
     all(rem(values, 1) == 0) && ...
     all(~isinf(values)) && ...
     all(values > 0) && ...
     (numel(values) == 1);
