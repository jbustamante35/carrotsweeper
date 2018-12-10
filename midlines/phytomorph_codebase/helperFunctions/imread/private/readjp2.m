function [A, map] = readjp2(filename, varargin)
%READJP2 Read image data from JPEG 2000 files.
%   A = READJP2(FILENAME) reads image data from a JPEG file.
%   A is a 2-D grayscale or 3-D RGB image whose type depends on the
%   bit-depth of the image (logical, uint8, uint16, int8, int16).
%
%   See also IMREAD.

%   Copyright 2008-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.11 $  $Date: 2011/05/17 02:28:03 $

% Error if anything other than a filename was passed.

% Check input arguments
options = parse_args(varargin{:});
if (~isNonnegativeFiniteIntegerScalar(options.reductionlevel))    
    error(message('MATLAB:imagesci:readjp2:badReductionLevel'))
end
options.pixelregion = process_region(options.pixelregion);
if ~isa(options.v79compatible, 'logical')
    error(message('MATLAB:imagesci:readjp2:badV79Compatible'))
end

% JPEG2000 is not supported on Solaris.
if (isequal(computer(), 'SOL64'))
    error(message('MATLAB:imagesci:readjp2:unsupportedPlatform'))
end

% Setup default options.
options.useResilientMode = false;  % default is fast mode

% Call the interface to the Kakadu library.
try
	A = readjp2c(filename,options);

catch firstException
	
	switch firstException.identifier
		case 'MATLAB:imagesci:jp2adapter:ephMarkerNotFollowingPacketHeader'

		    % Try resilient mode.  
			options.useResilientMode = true;
			try
				A = readjp2c(filename,options);

				% Ok we succeeded.  Issue a warning to the user that their
				% file might have some problems.  
				warning(message('MATLAB:imagesci:readjp2:ephMarkerNotFollowingPacketHeader', filename, firstException.message));

			catch secondException
				% Ok it's hopeless, just give up.
				rethrow(firstException);	
			end

		otherwise
			% We don't know what to try.  Give up.
			rethrow(firstException);	
	end


end
map = [];

function args = parse_args(varargin)
%PARSE_ARGS  Convert input arguments to structure of arguments.

args.reductionlevel = 0;
args.pixelregion = [];
args.v79compatible = false;

params = {'reductionlevel', 'pixelregion', 'v79compatible'};

p = 1;
while (p <= nargin)
    
    argp = varargin{p};
    if (ischar(argp))
        
        idx = find(strncmpi(argp, params, numel(argp)));
        
        if (isempty(idx))
            error(message('MATLAB:imagesci:readjp2:unknownParam', argp))
        elseif (numel(idx) > 1)
            error(message('MATLAB:imagesci:readjp2:ambiguousParam', argp))
        end
        
        if (p == nargin)
            error(message('MATLAB:imagesci:readjp2:missingValue', argp))
        end
        
        args.(params{idx}) = varargin{p + 1};
        p = p + 2;
        
    else
        
        error(message('MATLAB:imagesci:readjp2:paramType'))
        
    end
            
end


function region_struct = process_region(region_cell)
%PROCESS_PIXELREGION  Convert a cells of pixel region info to a struct.

region_struct = struct([]);
if isempty(region_cell)
    % Not specified in call to readjp2.
    return;
end

if ((~iscell(region_cell)) || (numel(region_cell) ~= 2))
    error(message('MATLAB:imagesci:readjp2:pixelRegionCell'))
end

for p = 1:numel(region_cell)
    
    checkIntegers(region_cell{p});
    
    if (numel(region_cell{p}) == 2)
        
        start = max(0, region_cell{p}(1) - 1);
        stop = region_cell{p}(2) - 1;
        
    else
        
        error(message('MATLAB:imagesci:readjp2:tooManyPixelRegionParts'));
        
    end
        
    if (start > stop)
        error(message('MATLAB:imagesci:readjp2:badPixelRegionStartStop'))
    end

    region_struct(p).start = start;
    region_struct(p).stop = stop;

end


function checkIntegers(inputValue)

if (~isnumeric(inputValue) || ...
    any((rem(inputValue, 1) ~= 0) & ~(isinf(inputValue))) || ...
    any(inputValue <= 0))
    
    error(message('MATLAB:imagesci:readjp2:regionPartNotNumeric'));

end


function tf = isNonnegativeFiniteIntegerScalar(values)

tf = isnumeric(values) && ...
     all(rem(values, 1) == 0) && ...
     all(~isinf(values)) && ...
     all(values >= 0) && ...
     (numel(values) == 1);

