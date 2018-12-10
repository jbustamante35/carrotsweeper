function writehdf(data, map, filename, varargin)
%WRITEHDF Write an raster-image HDF file to disk.
%   WRITEHDF(X,MAP,FILENAME) writes the indexed image X,MAP
%   to the file specified by the string FILENAME.
%
%   WRITEHDF(GRAY,[],FILENAME) writes the grayscale image GRAY
%   to the file.
%
%   WRITEHDF(RGB,[],FILENAME) writes the truecolor image
%   represented by the M-by-N-by-3 array RGB.
%
%   WRITEHDF(..., 'compression', COMP) uses the compression
%   type indicated by the string COMP.  COMP can be 'none',
%   'rle' (for grayscale or indexed images), or 'jpeg'
%   (for grayscale or rgb images).  The default value is
%   'none'.
%
%   WRITEHDF(..., 'compression', 'jpeg', 'quality', QUAL)
%   specifies the quality factor to use with JPEG
%   compression.  100 is best, 0 is worst, and 75 is
%   the default.  Lower values generally result in smaller
%   files.
%
%   WRITEHDF(..., 'writemode', MODE) either writes over
%   an existing HDF file (if MODE is 'overwrite') or appends
%   the image to the existing HDF file (if MODE is 'append').
%   The default is 'overwrite'.

%   Steven L. Eddins, August 1996
%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.7 $  $Date: 2011/05/17 02:28:14 $

if (ndims(data) > 3)
    error(message('MATLAB:imagesci:writehdf:tooManyDims', ndims( data )));
end

% 1 component or 3 components?
ncomp = size(data,3);
if ((ncomp ~= 1) && (ncomp ~= 3))
    error(message('MATLAB:imagesci:writehdf:wrongNumberOfComponents', ncomp));
end

% HDF library expects uint8 data
if (~isa(data,'double') && ~isa(data,'single') && ~isa(data,'uint8'))
    error(message('MATLAB:imagesci:writehdf:badImageClass'));
end

% Set default parameters
compressionType = 'none';
quality = 75;
baseline = 1;
writeMode = 'overwrite';

% Process param/value pairs
propStrings = {'compression', ...
               'quality    ', ...
               'baseline   ', ...
               'writemode  '};
    
for k = 1:2:length(varargin)
    prop = lower(varargin{k});
    if (~ischar(prop))
        error(message('MATLAB:imagesci:writehdf:badParameter'));
    end
    idx = find(strncmp(prop,propStrings,numel(prop)));
    if (isempty(idx))
        error(message('MATLAB:imagesci:writehdf:unrecognizedParameter', prop));
    elseif (length(idx) > 1)
        error(message('MATLAB:imagesci:writehdf:ambiguousParameter', prop));
    end
    
    prop = deblank(propStrings{idx});
    
    switch prop
    case 'compression'
        compressionType = varargin{k+1};
        if (~strcmp(compressionType, 'none') && ...
                    ~strcmp(compressionType, 'rle') && ...
                    ~strcmp(compressionType, 'jpeg'))
            error(message('MATLAB:imagesci:writehdf:badCompressionValue'));
        end
        
    case 'quality'
        quality = varargin{k+1};
        if ((quality < 0) || (quality > 100))
            error(message('MATLAB:imagesci:writehdf:badQualityValue'));
        end
        
    case 'baseline'
        baseline = varargin{k+1};
        if ((baseline ~= 0) && (baseline ~= 1))
            error(message('MATLAB:imagesci:writehdf:badBaselineValue'));
        end
        
    case 'writemode'
        writeMode = varargin{k+1};
        if (~strcmp(writeMode, 'append') && ...
                    ~strcmp(writeMode, 'overwrite'))
            error(message('MATLAB:imagesci:writehdf:badWriteModeValue'));
        end
        
    end
    
end
        
FAIL = -1;
if (ncomp == 1)
    hdf('DFR8', 'restart');
    
    % Convert data to uint8 if necessary
    if (~isa(data, 'uint8'))
        if (isempty(map))
            % Assume intensity image in range [0, 1].
            data = uint8(255*data);
        else
            % Assume indexed image
            data = uint8(data-1);
        end
    end

    if (~isempty(map) && ~isa(map, 'uint8'))
        map = uint8(255*map);
    end
    
    
    if (~isempty(map))
        if (size(map,1) < 256)
            map(256,3) = 0;
        elseif (size(map,1) > 256)
            error(message('MATLAB:imagesci:writehdf:tooManyColormapEntries'));
        end
        status = hdf('DFR8', 'setpalette', map');
        if (status == FAIL)
            error(message('MATLAB:imagesci:writehdf:libhdfError', hdferror));
        end
    else
        status = hdf('DFR8', 'setpalette', []);
        if (status == FAIL)
            error(message('MATLAB:imagesci:writehdf:libhdfError', hdferror));
        end
    end
    
    if (strcmp(compressionType, 'jpeg'))
        status = hdf('DFR8', 'setcompress', 'jpeg', quality, baseline);
        if (status == FAIL)
            error(message('MATLAB:imagesci:writehdf:libhdfError', hdferror));
        end
    else
        hdf('DFR8', 'setcompress', compressionType);
    end
    
    if (strcmp(writeMode, 'append'))
        status = hdf('DFR8', 'addimage', filename, data', compressionType);
        
    else
        status = hdf('DFR8', 'putimage', filename, data', compressionType);
    end
    
    if (status == FAIL)
        error(message('MATLAB:imagesci:writehdf:libhdfError', hdferror));
    end
    
else
    % 3 components

    if (strcmp(compressionType, 'rle'))
        error(message('MATLAB:imagesci:writehdf:rleNotSupportedForRgb'));
    end

    if ((isa(data,'double')) || (isa(data,'single')))
        data = uint8(data*255);
    end

    hdf('DF24', 'restart');

    data = permute(data, [3 2 1]);  % HDF API uses C-style element ordering
    
    if (~isempty(map))
        warning(message('MATLAB:imagesci:writehdf:colormapIgnoredForRgb'));
    end
    
    if (strcmp(compressionType, 'jpeg'))
        status = hdf('DF24', 'setcompress', 'jpeg', quality, baseline);
    else
        status = hdf('DF24', 'setcompress', compressionType);
    end
    if (status == FAIL)
        error(message('MATLAB:imagesci:writehdf:libhdfError', hdferror));
    end
            
    status = hdf('DF24', 'setil', 'pixel');
    if (status == FAIL)
        error(message('MATLAB:imagesci:writehdf:libhdfError', hdferror));
    end
    
    if (strcmp(writeMode, 'append'))
        status = hdf('DF24', 'addimage', filename, data);
        
    else
        status = hdf('DF24', 'putimage', filename, data);
    end
    
    if (status == FAIL)
        error(message('MATLAB:imagesci:writehdf:libhdfError', hdferror));
    end
    
end

%%%
%%% Function hdferror
%%%
function str = hdferror

str = hdf('HE', 'string', hdf('HE', 'value', 1));
