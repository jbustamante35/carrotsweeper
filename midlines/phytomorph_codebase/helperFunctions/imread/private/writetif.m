function writetif(data, map, filename, varargin)
%WRITETIF Write a TIF file to disk.
%   WRITETIF(X,MAP,FILENAME) writes the indexed image X,MAP
%   to the file specified by the string FILENAME.
%
%   WRITETIF(GRAY,[],FILENAME) writes the grayscale image GRAY
%   to the file.
%
%   WRITETIF(RGB,[],FILENAME) writes the truecolor image
%   represented by the M-by-N-by-3 array RGB.
%
%   WRITETIF(..., 'compression', COMP) uses the compression
%   type indicated by the string COMP.  COMP can be 'packbits',
%   'ccitt', 'fax3', 'fax4', or 'none'.  'ccitt', 'fax3', and 'fax4'
%   are allowed for logical inputs only.  'packbits' is the default.
%
%   WRITETIF(..., 'description', DES) writes the string contained
%   in DES to the TIFF file as an ImageDescription tag.
%
%   WRITETIF(..., 'resolution', XYRes) uses the scalar in XYRes
%   for the XResolution and YResolution tags.
%
%   WRITETIF(..., 'rowsperstrip', RPS) uses the scalar in RPS for
%   the RowsPerStrip tag.  Specifying too large a value will result in
%   RowsPerStrip being equal to the number of rows in the image.
%
%   WRITETIF(..., 'colorspace', CS) writes a TIFF file using the
%   specified colorspace, either 'rgb', 'icclab', or 'cielab'.  The input
%   image array must be M-by-N-by-3.
%
%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.11 $  $Date: 2011/05/17 02:28:21 $

if (ndims(data) > 3)
    error(message('MATLAB:imagesci:writetif:tooManyDims', ndims( data )));
end

% 1 component or 3 or 4 (cmyk) components ?
ncomp = size(data,3);
if ((ncomp ~= 1) && (ncomp ~= 3) && (ncomp ~= 4))
    error (message('MATLAB:imagesci:writetif:invalidComponentsNumber', ncomp));
end


[compression, description, resolution, rowsperstrip, writemode, colorspace] ...
    = parse_param_value_pairs ( data, map, varargin{:} );


% JPEG compression requires RowsPerStrip to be a multiple of 8.  It 
% also cannot use logical data.  The library catches this on most, 
% but not all platforms.
if strcmpi(compression,'jpeg')
    if (mod(rowsperstrip,8) ~= 0)
        error(message('MATLAB:imagesci:writetif:badRowsPerStripForJpegCompression'));
    end 
    if isa(data,'logical') 
        error(message('MATLAB:imagesci:writetif:badDatatypeForJpegCompression'));
    end 
end

if (ndims(data) == 3) && (size(data,3) == 3) && ...
        (isequal(colorspace, 'cielab') || isequal(colorspace, 'icclab'))
    if isa(data, 'double')
        % Convert to 8-bit ICCLAB.
        data(:,:,1) = round(data(:,:,1) * 255 / 100);
        data(:,:,2:3) = round(data(:,:,2:3) + 128);
        data = uint8(data);
    end
    
    if isequal(colorspace, 'cielab')
        % Convert to "munged" cielab values before writing them with
        % wtifc.
        data = icclab2cielab(data);
    end
end

% Logical data with more than 1 channel are not supported
if islogical(data) && size(data, 3) ~= 1
    error(message('MATLAB:imagesci:writetif:unsupportedData'));
end


if(isa(data,'double'))
    if (~isempty(map))
        % Zero based indexing for colormaps
        data = uint8(data - 1);
    else
        % Clip and scale double data
        data = min(1, max(0, data));
        data = uint8(255 * data);
    end
elseif ~(isa(data, 'logical') || isa(data, 'uint8') || isa(data, 'uint16'))
    error(message('MATLAB:imagesci:writetif:unsupportedDataType', class( data )));    
end


if (~isempty(map) && ~isa(map,'uint16'))
    map = uint16(65535 * map);
end

% Now that everything else is verified, check that it's possible to
% open a file here.  If it's not, fail.  If it is, this call to FOPEN
% will either create a new file or open an existing file.  Let the TIFF C
% code decide whether to truncate it.
fid = fopen(filename, 'a');
if (fid < 0)
    error(message('MATLAB:imagesci:writetif:fileOpen', filename));
else
    fclose(fid);
end

%
% Pack up the required tags into a single structure.
required_tags.compression     = lower(compression);
required_tags.description     = description;
required_tags.resolution      = resolution;
required_tags.rowsperstrip    = rowsperstrip;
required_tags.imageHeight     = size(data,1);
required_tags.imageWidth      = size(data,2);
required_tags.samplesPerPixel = size(data,3);

wtifc(data, map, filename, writemode, colorspace, required_tags);



%===============================================================================
function [comp, desc, res, rps, wmode, colorspace] = parse_param_value_pairs ( data, map, varargin )
%
% comp = compression
% desc = description
% res = resolution
% rps = rowsperstrip
% wmode = writemode

desc = '';
comp= 'packbits';
res = 72;
rps = -1;
wmode = 1;
colorspace = 'rgb';


if (islogical(data) && (ndims(data) == 2) && isempty(map))
    comp= 'ccitt';
end

% Process param/value pairs
paramStrings = {'compression  ', ...
                'description  ', ...
                'resolution   ', ...
                'writemode    ', ...
                'colorspace   ', ...
                'rowsperstrip ' };
    
for k = 1:2:length(varargin)
    param = lower(varargin{k});
    if (~ischar(param))
        error(message('MATLAB:imagesci:writetif:badParameterName'));
    end
    idx = find(strncmp(param, paramStrings, numel(param)));
    if (isempty(idx))
        error(message('MATLAB:imagesci:writetif:unrecognizedParameter', param));
    elseif (length(idx) > 1)
        error(message('MATLAB:imagesci:writetif:ambiguousParameter', param));
    end
    
    param = deblank(paramStrings{idx});
    
    switch param
    case 'compression'
        comp = varargin{k+1};
        if (~ischar(comp))
            error(message('MATLAB:imagesci:writetif:badCompression'));
        end
        
    case 'description'
        desc = varargin{k+1};
        if (~ischar(desc))
            error(message('MATLAB:imagesci:writetif:badDescription'));
        end
        
    case 'resolution'
        res = varargin{k+1};
        
    case 'rowsperstrip'
        rps = varargin{k+1};
        if ~isnumeric(rps)
            error(message('MATLAB:imagesci:writetif:badRowsPerStrip'));
        end
        if isempty(rps)
            error(message('MATLAB:imagesci:writetif:emptyRowsPerStrip'));
        end
        if rps < 1
            error(message('MATLAB:imagesci:writetif:rowsPerStripTooSmall'));
        end
        
    case 'writemode'
        switch lower(varargin{k+1})
            case 'overwrite'
                wmode = 1;
                
            case 'append'
                wmode = 0;
                
            otherwise
                error(message('MATLAB:imagesci:writetif:badWriteMode'));
        end
        
    case 'colorspace'
        colorspace = varargin{k+1};
        if ~ischar(colorspace)
            error(message('MATLAB:imagesci:writetif:colorspaceInvalidClass'));
        end
        colorspace = lower(colorspace);
        if ~ismember(colorspace, {'rgb', 'cielab', 'icclab'})
            error(message('MATLAB:imagesci:writetif:colorspaceUnrecognizedString'));
        end

        if (ndims(data) ~= 3) && (size(data, 3) ~= 3)
            warning(message('MATLAB:imagesci:writetif:ignoredColorspaceInput'));
        end

    end
end

