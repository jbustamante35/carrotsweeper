function writepng(data, map, filename, varargin)
%WRITEPNG Write a PNG file to disk.
%   WRITEPNG(I,[],FILENAME) writes the grayscale image I
%   to the file specified by the string FILENAME.
%
%   WRITEPNG(RGB,[],FILENAME) writes the truecolor image
%   represented by the M-by-N-by-3 array RGB.
%
%   WRITEPNG(X,MAP,FILENAME) writes the indexed image X with
%   colormap MAP.  The resulting file will contain the equivalent
%   truecolor image.
%
%   WRITEPNG(...,PARAM,VAL,...) sets the specified parameters.
%
%   See also IMREAD, IMWRITE, IMFINFO.

%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.1.6.9 $  $Date: 2011/05/17 02:28:18 $

if (ndims(data) > 3)
    error(message('MATLAB:imagesci:writepng:tooManyDims', ndims( data )));
end

% Color type values (as in PNG library defs)
PNG_COLOR_TYPE_GRAY = 0;
PNG_COLOR_TYPE_RGB = 2;
PNG_COLOR_TYPE_PALETTE = 3;
PNG_COLOR_TYPE_GRAY_ALPHA = 4;
PNG_COLOR_TYPE_RGB_ALPHA = 6;

% Set default parameters
bitdepth = [];
sigbits = [];
interlace = 'none';
transparency = [];
alpha = [];
background = [];
gamma = [];
chromaticities = [];
xres = [];
yres = [];
resunit = [];
textchunks = cell(0,2);
imagemodtime = [];

% Process param/value pairs
propStrings = {'interlacetype  ', ...
    'transparency   ', ...
    'bitdepth       ', ...
    'significantbits', ...
    'alpha          ', ...
    'background     ', ...
    'gamma          ', ...
    'chromaticities ', ...
    'xresolution    ', ...
    'yresolution    ', ...
    'resolutionunit ', ...
    'title          ', ...
    'author         ', ...
    'description    ', ...
    'copyright      ', ...
    'creationtime   ', ...
    'software       ', ...
    'disclaimer     ', ...
    'warning        ', ...
    'source         ', ...
    'comment        ', ...
    'imagemodtime   '};

for k = 1:2:length(varargin)
    prop = lower(varargin{k});
    if (~ischar(prop))
        error(message('MATLAB:imagesci:writepng:parameterNotString'));
    end
    idx = find(strncmp(prop, propStrings, numel(prop)));
    if (isempty(idx))
        keyword = varargin{k};
        textItem = varargin{k+1};
        keyword = CheckKeyword(keyword);
        textItem = CheckTextItem(textItem);
        textchunks{end+1,1} = keyword;
        textchunks{end,2} = textItem;
    
    elseif (length(idx) > 1)
        error(message('MATLAB:imagesci:writepng:ambiguousParameter', prop));
        
    else
        prop = deblank(propStrings{idx});
        switch prop
        case 'bitdepth'
            bitdepth = varargin{k+1};
            
        case 'significantbits'
            sigbits = varargin{k+1};
            
        case 'interlacetype'
            interlace = varargin{k+1};
            
        case 'transparency'
            transparency = varargin{k+1};
            
        case 'alpha'
            alpha = varargin{k+1};
            
        case 'background'
            background = varargin{k+1};
            
        case 'gamma'
            gamma = varargin{k+1};
            
        case 'chromaticities'
            chromaticities = varargin{k+1};
            
        case 'xresolution'
            xres = varargin{k+1};
            
        case 'yresolution'
            yres = varargin{k+1};
            
        case 'resolutionunit'
            resunit = varargin{k+1};
            
        case 'title'
            keyword = 'Title';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'author'
            keyword = 'Author';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'description'
            keyword = 'Description';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'copyright'
            keyword = 'Copyright';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'creationtime'
            keyword = 'Creation Time';
            if (ischar(varargin{k+1}))
                textItem = datestr(datenum(varargin{k+1}), 0);
            else
                textItem = datestr(varargin{k+1}, 0);
            end
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'software'
            keyword = 'Software';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'disclaimer'
            keyword = 'Disclaimer';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'warning'
            keyword = 'Warning';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'source'
            keyword = 'Source';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'comment'
            keyword = 'Comment';
            textItem = CheckTextItem(varargin{k+1});
            textchunks{end+1,1} = keyword;
            textchunks{end,2} = textItem;
            
        case 'imagemodtime'
            try
                imagemodtime = fix(datevec(varargin{k+1}));
            catch
                error(message('MATLAB:imagesci:writepng:invalidImageModTime'));
            end
            
            if (numel(imagemodtime) > 6)
                error(message('MATLAB:imagesci:writepng:tooMuchImageModTimeData'))
            end
            
        end
    end
    
end

if ((ndims(data) > 3) || (~ismember(size(data,3), [1 3])))
    error(message('MATLAB:imagesci:writepng:wrongImageDimensions'));
end

if (~ismember({class(data)}, {'double', 'single', 'logical', 'uint8', 'uint16'}))
    error(message('MATLAB:imagesci:writepng:unsupportedImageClass'));
end

if (~isempty(alpha) && ((size(alpha,1) ~= size(data,1)) || ...
                    (size(alpha,2) ~= size(data,2))))
    error(message('MATLAB:imagesci:writepng:alphaDoesNotMatchImage'));
end

%
% Identify color type
%
isTruecolor = (size(data,3) == 3);
paletteUsed = ~isempty(map) && ~isTruecolor;
colorUsed = paletteUsed || isTruecolor;
alphaUsed = ~isempty(alpha);
colortype = paletteUsed + 2*colorUsed + 4*alphaUsed;
if (colortype == 7)
    error(message('MATLAB:imagesci:writepng:alphaNotSupportedForIndexed'));
end

%
% Set default bitdepth if not specified
%
if (isempty(bitdepth))
    switch class(data)
        case 'logical'
            bitdepth = 1;

        case {'uint8', 'double', 'single'}
            bitdepth = 8;

        case 'uint16'
            bitdepth = 16;
    end
end

%
% Validate bitdepth
%
switch colortype
    case PNG_COLOR_TYPE_GRAY
        if (~ismember(bitdepth, [1 2 4 8 16]))
            error(message('MATLAB:imagesci:writepng:invalidGrayscaleBitDepth'));
        end
        
    case PNG_COLOR_TYPE_RGB
        if (~ismember(bitdepth, [8 16]))
            error(message('MATLAB:imagesci:writepng:invalidRgbBitDepth'));
        end
        
    case PNG_COLOR_TYPE_PALETTE
        if (~ismember(bitdepth, [1 2 4 8]))
            error(message('MATLAB:imagesci:writepng:invalidIndexedBitDepth'));
        end
        
    case PNG_COLOR_TYPE_GRAY_ALPHA
        if (~ismember(bitdepth, [8 16]))
            error(message('MATLAB:imagesci:writepng:invalidGrayscaleAlphaBitDepth'));
        end
        
    case PNG_COLOR_TYPE_RGB_ALPHA
        if (~ismember(bitdepth, [8 16]))
            error(message('MATLAB:imagesci:writepng:invalidRgbAlphaBitDepth'));
        end
end

%
% Scale image if necessary to match requested bitdepth
%
switch class(data)
    case {'double', 'single'}
        if (colortype == PNG_COLOR_TYPE_PALETTE)
            data = data - 1;
            data = uint8(data);
        
        else
            % Grayscale or RGB; clamp data to [0,1] dynamic range before
            % scaling, rounding, and casting.
            data = max(min(data,1),0);
            switch bitdepth
                case 8
                    data = uint8(255*data);
                    
                case 16
                    data = uint16(65535*data);
                    
                case 4
                    data = uint8(15*data);
                    
                case 2
                    data = uint8(3*data);
                    
                case 1
                    data = uint8(data ~= 0);
            end
        end
        
    case 'uint8'
        if (colortype == PNG_COLOR_TYPE_PALETTE)
            % Nothing to do
            
        else
            switch bitdepth
                case 16
                    data = uint16(data);
                    data = bitor(bitshift(data,8),data);
                    
                case 8
                    % Nothing to do
                    
                case 4
                    data = bitshift(data,-4);
                    
                case 2
                    data = bitshift(data,-6);
                    
                case 1
                    % Nothing to do
            end
        end
        
    case 'uint16'
        switch bitdepth
            case 16
                % Nothing to do
                
            case 8
                data = uint8(bitshift(data,-8));
                
            case 4
                data = uint8(bitshift(data,-12));
                    
            case 2
                data = uint8(bitshift(data,-14));
                    
            case 1
                data = uint8(data ~= 0);
        end
end

if (ismember(colortype, [PNG_COLOR_TYPE_GRAY_ALPHA, ...
                        PNG_COLOR_TYPE_RGB_ALPHA]))
    %
    % Scale alpha data if necessary to match data class
    %
    switch bitdepth
        case 8
            switch class(alpha)
                case {'double', 'single'}
                    alpha = max(min(alpha,1),0);
                    alpha = uint8(255 * alpha);
                    
                case 'uint16'
                    alpha = uint8(bitshift(alpha, -8));
                    
                case 'uint8'
                    % nothing to do
                    
                otherwise
                    error(message('MATLAB:imagesci:writepng:bad8bitAlphaClass'));
            end
            
        case 16
            switch class(alpha)
                case {'double', 'single'}
                    alpha = max(min(alpha,1),0);
                    alpha = uint16(65535 * alpha);
                    
                case 'uint16'
                    % nothing to do
                    
                case 'uint8'
                    alpha = uint16(alpha);
                    alpha = bitor(bitshift(alpha, 8), alpha);
                    
                otherwise
                    error(message('MATLAB:imagesci:writepng:badAlphaClass'));
            end
    end
end

% Be friendly about specifying resolutions
if (~isempty(xres) && isempty(yres))
    yres = xres;

elseif (~isempty(yres) && isempty(xres))
    xres = yres;
end

if (~isempty(xres) && isempty(resunit))
    resunit = 'unknown';
end

if (isempty(xres) && isempty(yres) && ~isempty(resunit))
    error(message('MATLAB:imagesci:writepng:resolutionsRequired'));
end
        
pngwritec(data, map, filename, colortype, bitdepth, ...
                sigbits, alpha, interlace, ...
                transparency, background, gamma, ...
                chromaticities, xres, yres, ... 
                resunit, textchunks, imagemodtime);


function out = CheckKeyword(in)
%CheckKeyword
%   out = CheckKeyWord(in) checks the validity of the input text chunk keyword.

if (isempty(in))
    error(message('MATLAB:imagesci:writepng:emptyTextChunkKeyword'))
end
if ((in(1) == 32) || (in(end) == 32))
    error(message('MATLAB:imagesci:writepng:paddedTextChunkKeyword'));
end
if (numel(in) > 80)
    error(message('MATLAB:imagesci:writepng:tooMuchKeywordData'));
end
if (any(~ismember(in,[32:126 161:255])))
    error(message('MATLAB:imagesci:writepng:invalidCharsInTextChunkKeyword'));
end

out = in;


function out = CheckTextItem(in)
%CheckTextItem
%   out = CheckTextItem(in) strips out control characters from text; PNG spec
%   discourages them.  It also replaces [13 10] by 10; then it replaces 13 
%   by 10.  The PNG spec says newlines must be represented by a single 10.

if (~ischar(in))
    error(message('MATLAB:imagesci:writepng:invalidTextChunk'));
end

out = in;
out = strrep(out, char([13 10]), char(10));
out = strrep(out, char(13), char(10));
badChars = find((out < 32) & (out ~= 10));
if (~isempty(badChars))
    warning(message('MATLAB:imagesci:writepng:changedTextChunk'));
    out(badChars) = [];
end
